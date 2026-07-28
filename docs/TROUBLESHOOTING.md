# Troubleshooting log

Real debugging episodes, their root cause, and the fix — kept as an engineering
record.

---

## 2026-07-28 — Per-module AI (`ai_interpret_modules`) 500s: Ollama `num_ctx` reload thrash

### Symptom

On a full run with the local 30B model (`huihui_ai/qwen3-abliterated:30b-a3b`,
~25 GB unified-memory Mac), the per-category AI, executive summary, and
cross-category synthesis all completed — but **every one of the 9 per-module
AI interpretations failed instantly** with:

```
WARNING: module AI for holistic-synthesis failed: Ollama error: 500 Server Error:
  Internal Server Error for url: http://localhost:11434/api/chat
WARNING: module AI for immunogenetics failed: Ollama error: 500 ...
... (all 9 modules) ...
Per-module AI: 0 module interpretations generated
```

The report still generated correctly (per-module AI failures are non-fatal),
but "AI on all tiers" produced nothing.

### Diagnosis (reproduced, not guessed)

A controlled reproduction against the same model showed the failure was tied to
**`num_ctx`**, not the prompt content:

| call | `num_ctx` | `think` | result |
|------|-----------|---------|--------|
| 1    | 16384     | False   | ✅ OK (9s) |
| 2    | 4096      | False   | ❌ 500 (6s) |
| 3    | 4096      | False   | ❌ 500 (0s) |
| 4    | 16384     | False   | ✅ OK (7s) |

Only the context size that the model was **already loaded with** succeeded;
requesting a *different* `num_ctx` failed. Root cause:

- `num_ctx` is a **load-time** parameter in Ollama. Changing it between requests
  forces the model to be **reloaded**.
- Reloading a large (30B) model on a 25 GB unified-memory machine intermittently
  **OOMs during the reload**, which surfaces as an HTTP **500**.

In the real run the sequence was: per-category loaded the model at
`num_ctx=4096` → executive summary + cross-category reloaded it to `16384`
(succeeded) → the per-module calls requested `4096` again → **reload → 500**,
and every subsequent call 500'd because the runner was in a bad state.

Why per-category survived but per-module didn't: per-category uses the default
(`think=None`) and ran *first* (loading the model), whereas the per-module calls
requested a smaller context *after* the 16384 exec/cross calls — the size change
was the trigger.

### Fix

1. **One context size for every Ollama call.** Introduced `AI_NUM_CTX = 16384`
   in `analyze.py` and routed per-category, executive summary, cross-category,
   and per-module calls all through it, so the model **loads once and never
   reloads**. 16384 comfortably fits the largest prompt (cross-category, capped
   at 48k chars ≈ 12k tokens + output).
2. **`keep_alive: "15m"`** on the request so the model stays resident across the
   many sequential tier-2 calls.
3. **Retry on transient 5xx** (2 attempts, short backoff) in `call_ollama` as a
   safety net, since a single momentary reload/OOM should not blank a section.

### Verification

Three consecutive `think=False` calls at the standardized context now succeed
(9s, then 1s, 1s — the model stays resident, no reload):

```
AI_NUM_CTX = 16384
  call 0 (default num_ctx, think=False): OK 9s
  call 1 (default num_ctx, think=False): OK 1s
  call 2 (default num_ctx, think=False): OK 1s
```

### Related earlier fix (same file)

`think=False` was itself added earlier to fix a *different* failure: reasoning
models (qwen3 / abliterated variants) emit `<think>…</think>` blocks that
`call_ollama` strips; with a small `num_predict` the whole budget could be
consumed inside `<think>`, stripping to an **empty** section. The synthesis /
interpretation calls disable reasoning so the budget produces the visible
answer, plus a salvage guard so stripping never yields `""`.
