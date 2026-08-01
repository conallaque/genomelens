# GenomeLens — Full Methodology

Complete mathematical derivation, citations, and worked examples behind every
number GenomeLens reports. See the main README for a plain-language overview;
this document is the technical appendix for readers who want to verify or extend
the model.

---

### The math — with a plain-English translation under every equation

Written like a methods section, but **each equation has a translation underneath** so
you don't need a stats background to follow what's happening. Equations use GitHub's
math blocks; every symbol is also spelled out in words.

**Notation.** For a finding `i`: `λ` = how much one year of healthy life is worth in
dollars; `p_i` = the chance the condition actually happens; `RRR_i` = how much acting
cuts that risk; `COI_i` = the lifetime cost of the illness; `q_i` = healthy life-years
(QALYs) preserved by acting; `h_i` = a confidence discount (1 for confirmed findings,
less for computational predictions); `T_i` = the years it plays out over; `r` = the
discount rate; `C_test` = the cost of the DNA test.

**1 · Discounting — value the future a little less than today**

```math
D(T,r)=\frac{1}{T}\left(\frac{1}{2}+\sum_{t=1}^{T}\frac{1}{(1+r)^{t}}\right),\qquad r=0.03
```

*In plain English:* a dollar — or a healthy year — decades from now is worth less than
one today. This shrinks future costs **and** future health benefits by 3% per year so we
never overstate things that pay off far away. It's the standard, conservative convention.

**2 · What acting on a finding is expected to be worth**

```math
\Delta\text{Cost}_i = p_i\cdot\text{RRR}_i\cdot\text{COI}_i\cdot D(T_i,r)\cdot h_i
```
```math
\Delta\text{QALY}_i = p_i\cdot\text{RRR}_i\cdot q_i\cdot D(T_i,r)\cdot h_i
```

*In plain English:* the benefit of acting = **how likely the problem is** × **how much
acting reduces it** × **how costly/harmful it is** — discounted for time, and dialed back
if the finding is only a prediction. One line is money saved; the other is healthy years
gained.

**2b · Drug-response (pharmacogenomic) findings**

```math
\Delta\text{Cost}^{\text{PGx}}_i = P(\text{Rx})\cdot P(\text{ADR}\mid g)\cdot\text{RRR}\cdot C_\text{ADR}
```

*In plain English:* for a medication gene, the value = **chance you'll be prescribed the
drug** × **chance it harms you given your genotype** × **how much genotype-guided dosing
avoids that** × **cost of the bad reaction**. In short: the harm and money avoided by not
getting the wrong prescription.

**3 · Net Monetary Benefit — the bottom line for each finding**

```math
\text{NMB}_i = \lambda\cdot\Delta\text{QALY}_i + \Delta\text{Cost}_i - c^{\text{int}}_i
```

*In plain English:* turn the healthy years gained into dollars (at `λ`), add the money
saved, then subtract what acting costs (the screening, drug, or visit, `c_int`).
**Positive means it's worth doing.** We lead with this instead of the usual ratio because
it still behaves sensibly when a step actually *saves* money.

**4 · Value of Information — what the whole test is worth**

```math
\text{VOI}_{\text{ex-post}}=\sum_{i}\text{NMB}_i - C_\text{test}
```
```math
\text{VOI}_{\text{ex-ante}}=\sum_{k}\pi_k\,\text{NMB}_k - C_\text{test}
```
```math
\text{VOI}_{\Delta}=\sum_{i\,\in\,\text{WGS-only}}\text{NMB}_i
```

*In plain English:* add up the value of every finding, then subtract what the test cost.
**Ex-post** = what *your* results are worth. **Ex-ante** = what testing is worth to a
typical person *before* they know their results (`π_k` = how common finding *k* is in the
population). The third line (**VOI-Δ**) isolates how much *more* a full genome gives you
than a cheap chip — i.e. "is upgrading worth it?"

**5 · Honest uncertainty — the simulation behind the range**

```math
P_{\text{CE}}(\lambda)=\frac{1}{N}\sum_{j=1}^{N}\mathbb{1}\!\left[\text{NMB}^{(j)}(\lambda)>0\right],\qquad N=10{,}000
```

*In plain English:* every input above is uncertain, so instead of one number we run the
whole calculation **10,000 times**, each time drawing plausible values (chances from a
Beta distribution, costs from a Gamma distribution). That yields a **range** (a 95%
confidence interval) rather than false precision. The curve this produces — the CEAC —
reads as: *"at a given value of a healthy year, what's the probability this is actually
worth it?"*

**6 · Left-tail risk — how bad could this reasonably go?**

```math
\mathrm{VaR}_{95}=F^{-1}(0.05),\qquad
\mathrm{CVaR}_{95}=\mathbb{E}\!\left[V \mid V\le \mathrm{VaR}_{95}\right]
```

*In plain English:* an average hides the downside. **VaR** is the 5th-percentile outcome —
a bad-but-plausible case. **CVaR** (expected shortfall) is the *average of everything at or
below that*, i.e. "if it does go badly, how badly on average?" CVaR is the measure risk
managers prefer because, unlike VaR, it's *coherent* — diversifying can never make it look
artificially worse. This is the metric behind the tool's "left-tail risk" framing.

**7 · Health as depreciating capital (Grossman 1972)**

```math
H_{t+1}=H_t\bigl(1-\delta(a)\bigr)+I\cdot(1+\varepsilon),\qquad
\delta(a)=\delta_0 e^{g(a-20)}
```

*In plain English:* your health is a **stock**, like a machine or a portfolio. It wears out
a bit each year — and the wear rate `δ` accelerates as you age — while your effort `I`
(exercise, screening, treatment) tops it back up. Genomic information doesn't add health
directly; it makes each unit of effort **more efficient** (`ε`) because you target what
actually matters for you. The key consequence: an efficiency gain compounds over the years
you have left, so **the same information is worth more the younger you are** — the model
reproduces exactly that.

**8 · When to test — the option to wait (Dixit & Pindyck 1994)**

```math
V(T)=\frac{\text{VOI}\cdot\left(1-\tfrac{T}{H}\right)-C\,(1-c)^{T}}{(1+r)^{T}}
```

*In plain English:* sequencing gets cheaper every year, so why not wait? Because two things
cut the other way. First, **you can re-analyse stored data for free forever** — so improving
science reaches early testers too, and isn't a reason to delay. Second, **every year you
wait is a year you can't act** on what you'd have found (the `1 − T/H` term). For a genome
with real findings, waiting destroys more value than the price drop recovers → *test now*.
The model still says "wait" when the value is low relative to price — which is the honest
answer in that case.

**9 · The ceiling on information value — EVPI**

```math
\text{EVPI}=\mathbb{E}_{\theta}\!\left[\max_{a} \text{NB}(a,\theta)\right]-\max_{a}\mathbb{E}_{\theta}\!\left[\text{NB}(a,\theta)\right]\;\ge\;0
```

*In plain English:* imagine a clairvoyant who removes **all** remaining uncertainty. EVPI is
the most that clairvoyance could possibly be worth — therefore an upper bound on the value of
*any* further study or confirmatory test. Crucially, **a small EVPI is good news**: it means
the recommended actions stay optimal across nearly the whole uncertainty range, so the
decision is robust and more data wouldn't change it. Low EVPI means "act," not "weak analysis."

**10 · Risk preferences — why information is worth more than its average (Arrow–Pratt)**

```math
u(\text{CE})=\mathbb{E}\bigl[u(W+X)\bigr]
\;\;\Longrightarrow\;\;
\text{CE}\approx\mu-\tfrac{1}{2}\,\gamma\,\frac{\sigma^{2}}{W+\mu}
```

*In plain English:* most people would take a guaranteed \$900 over a coin-flip for \$2,000 —
that gap is **risk aversion** (`γ`). The **certainty equivalent** is the guaranteed amount
you'd accept instead of the uncertain outcome; the difference from the average is the **risk
premium** you'd pay to avoid uncertainty. Since genomic information *reduces* uncertainty
about your health, a risk-averse person values it **above** its expected dollar value — the
same logic that makes insurance markets exist (Arrow 1963).

**11 · Correcting the genetics before the economics**

```math
\text{odds}_{\text{pop}}=\frac{\text{odds}_{\text{lit}}}{\kappa},\qquad
\text{odds}_{\text{post}}=\text{odds}_{\text{pop}}\times \text{BF}_{\text{FH}},\qquad
\hat\beta_{\text{shrunk}}=\hat\beta\cdot\phi(|z|)
```

*In plain English:* two well-known statistical traps would silently inflate every dollar
figure, so both are corrected **before** the economics runs:

- **Ascertainment bias.** Classic penetrance numbers come from families studied *because*
  they had lots of cancer. Applying those to someone who found a variant incidentally
  overstates their risk — so the estimate is shrunk toward the population rate (`κ`), then
  updated on family history (`BF`).
- **Winner's curse.** A gene effect discovered *because* it cleared a significance cutoff is
  biased upward — the discovery sample's noise had to help it get found. Effect sizes are
  shrunk (`φ`) accordingly.

Both make the final numbers **smaller and more defensible** — which is the point.

**12 · From percentile to actual risk — the liability-threshold model (Falconer 1965)**

```math
T=\Phi^{-1}(1-K),\qquad
P(\text{affected}\mid \text{PRS})=1-\Phi\!\left(\frac{T-z_q\sqrt{R^2}}{\sqrt{1-R^2}}\right)
```

*In plain English:* "you're in the 95th percentile" is not a risk. This is the standard
way to turn it into one. Picture a hidden "liability" scale everyone sits on; you get the
disease if you cross a threshold `T` set by how common it is (`K`). Your score nudges you
along that scale, and the answer is how much of your remaining bell curve sits past the
line. **Crucially it returns absolute risk, not just relative risk** — a "3× higher risk"
of something rare is still rare, and that's the single most common way consumer genomics
misleads people.

*(A subtlety worth knowing: risk is convex here, so the median-PRS person sits slightly
below the population average while the average across everyone comes out exactly at `K`.
Verified numerically — see the docstring.)*

**13 · Penetrance is a curve, not a number — with competing risks (Fine & Gray 1999)**

```math
\text{CIF}(t)=\int_0^{t} S_{\text{all}}(u)\,h_{\text{disease}}(u)\,du
\;<\;1-e^{-\int_0^t h_{\text{disease}}}
```

*In plain English:* "55% lifetime risk" hides *when*, and ignores that **you can only die
once**. If something else takes you first, you never get the disease. The naive formula
(right side) quietly assumes nobody dies of anything else and therefore **overstates
risk** — by about 9 percentage points in our worked case. The left side is the honest
version: risk accumulated only over the years you're actually alive to accumulate it.

**14 · Longer lifespans make this worth *more*, not less**

Rising life expectancy pushes in the same direction twice: less competing mortality means
more genetic risk is actually realised, **and** a longer horizon means prevention compounds
for more years. Modelled explicitly:

| Scenario | Life expectancy | Realised risk | vs today |
|---|---|---|---|
| 2025 baseline | ~79 | 73.2% | — |
| 2050 projection | ~85 | 77.0% | +5.3% |
| Longevity-advance | ~95 | 81.0% | +10.8% |

*In plain English:* the model's default assumes today's lifespans. If medicine and applied
AI extend healthy life within our lifetimes — plausible — then the value reported here is an
**under**-estimate, not an over-estimate. The conservative choice is genuinely conservative.

**15 · Shrinking many noisy estimates at once (James–Stein / empirical Bayes)**

```math
\tau^{2}=\max\!\left(0,\ \widehat{\mathrm{Var}}(\hat\beta)-\overline{SE^{2}}\right),\quad
w_i=\frac{\tau^{2}}{\tau^{2}+SE_i^{2}},\quad
\beta_i^{*}=w_i\hat\beta_i+(1-w_i)\bar\beta
```

*In plain English:* when you hold many noisy estimates, taking each at face value is
provably worse than pulling them all toward the average — the famous James–Stein result.
The weight `w` is each estimate's reliability, so **the noisiest numbers get pulled the
hardest**. Same instinct as not trading on a single noisy signal.

**16 · Polygenic scores don't transfer across ancestries — so the model attenuates them**

```math
R^{2}_{\text{eff}}=R^{2}_{\text{EUR}}\times \rho_{\text{ancestry}},\qquad
\rho \approx 1.00\,/\,0.65\,/\,0.50\,/\,0.25
```

*In plain English:* most genetic studies were done in European-ancestry cohorts, so those
scores predict less well elsewhere — retaining roughly a quarter of their accuracy in
African-ancestry individuals. Because risk runs through the liability model, **a less
informative score automatically pulls your estimate back toward the population average**,
which is the statistically correct behaviour. This is a data-equity problem in genomics,
and the model surfaces it rather than burying it in a footnote.

**17 · Which uncertainty is worth resolving (EVPPI), and why people still don't test**

```math
\text{EVPPI}(\varphi)=\mathbb{E}_{\varphi}\!\left[\max_a \mathbb{E}_{\theta\mid\varphi}\text{NB}\right]-\max_a\mathbb{E}_{\theta}\text{NB},
\qquad 0\le\text{EVPPI}\le\text{EVPI}
```

*In plain English:* EVPPI asks, one input at a time, *"if I could know this perfectly,
would I actually decide differently?"* It tells you where a confirmatory test is worth
buying — and where more precision would be wasted money.

And a behavioural coda: under **prospect theory**, a certain up-front cost hurts about
**2.25×** more than an equivalent uncertain future gain feels good, while **hyperbolic
discounting** further shrinks benefits that arrive decades out. Together they explain the
gap between a test that is clearly worth it on paper and one people actually buy — an
*adoption* problem, not a valuation problem.

**18 · Asymmetric information — what the model deliberately does not price**

Genomic results create textbook **asymmetric information**: if you know your risk and an
insurer doesn't, that's **adverse selection** (Akerlof 1970; Rothschild–Stiglitz 1976); if the
insurer can price on it, fear of **genetic discrimination** deters testing, and socially
valuable information goes unacquired. In the US, **GINA (2008) covers health insurance and
employment but *not* life, disability, or long-term-care insurance.**

**19 · The welfare case for local analysis — and the number that decides it**

This is the economic argument for the tool's whole design, stated so it can be *falsified*
rather than asserted. The naive version is circular:

```math
S_{\text{local}} = B - 0 \;>\; S_{\text{central}} = B - C_{\text{privacy}}
```

That holds for any positive privacy cost **by construction**, because it assumes both
options deliver the same gross benefit. So the model instead **concedes that centralised
platforms are better at the analysis** (bigger reference panels, curated pipelines, expert
review) by a capability premium `π`:

```math
S_{\text{local}} = B_L - C_{\text{test}},\qquad
S_{\text{central}} = \underbrace{B_L(1+\pi)}_{\text{better analysis}} - C_{\text{test}} - \mathbb{E}[L_{\text{privacy}}]
```

```math
\textbf{Local wins} \iff \mathbb{E}[L_{\text{privacy}}] > \pi B_L
```

*In plain English:* going local costs you something real — a platform with more data may
genuinely analyse it better. Local only wins if the expected privacy cost is **bigger than
that sacrifice**. Note what cancels: `B_L` and the sequencing cost drop out of the
comparison entirely, because you pay the test either way.

**Exposure is an absorbing state.** A genome cannot be re-keyed after disclosure, so the
hazard applies every year the data sits in someone else's system, and one breach is enough
— forever:

```math
\mathbb{E}[L_{\text{privacy}}] = \underbrace{\bigl[1-(1-p)^{T}\bigr]}_{\text{at least one breach in }T\text{ yrs}} \cdot L \cdot (1+r)^{-T/2}
```

*(Not `p × T` — that double-counts years after the first breach, and exceeds 1 when
`pT > 1`. It's a survival problem, not a counting problem.)*

**The break-even.** Setting the two surpluses equal and solving gives the annual breach
probability at which a rational person is exactly indifferent:

```math
p^{*} = 1-\left(1-\frac{\pi B_L}{L\,(1+r)^{-T/2}}\right)^{1/T}
```

With the shipped defaults (`B_L` = \$25k, `π` = 15%, `L` = \$25k, `T` = 40, `r` = 3%):
**`p* ≈ 0.79% per year.`** Below that, centralised analysis genuinely wins; above it, local
does. **The model can and does conclude "centralised wins" at low breach risk** — which is
precisely what makes the opposite conclusion meaningful.

**Then check it against reality.** In 2023 a single credential-stuffing incident exposed
roughly **6.9 million of ~14 million** 23andMe users — about **49% in one event**, some
**60×** the 0.79% annual threshold. The empirical hazard clears the bar decisively.

**The bigger channel is access, not exposure — and this one is empirically grounded.**
Disclosure risk doesn't just impose a cost on people who upload; it stops people testing at
all, and a non-tester forfeits *100%* of the value. Weighting surplus by participation (`θ`):

```math
W = \theta \cdot S,\qquad
\underbrace{(\theta_L-\theta_C)\,S_{\text{local}}}_{\text{access channel}}
```

The participation gap is **not an assumption**. Miller & Tucker (2018, *Management Science*)
exploit state-level variation in US genetic-privacy law and find that regimes granting
patients **control** over their genetic data raise testing incidence by **+83%**, while
regimes that merely notify people of privacy risk and ask them to consent — *without*
granting control — **lower testing by 69%**. That contrast maps almost directly onto the
choice modelled here: local analysis is the maximal-control regime; uploading under a
terms-of-service click is the notice-without-control regime. Survey evidence agrees on
direction — NORC finds **~80%** of Americans hold privacy concerns about DNA testing, **~17%**
of non-testers name privacy as their reason, and **four in five** non-testers say they would
be more willing if privacy were assured.

The defaults used here (`θ_L` = 0.85 vs `θ_C` = 0.60, a 1.42× ratio) are **deliberately
more conservative than that literature**, whose implied control-vs-notice ratio is far
larger. Even so, the access channel is worth **\$6,175** — *larger than the entire
per-person privacy effect (\$3,923)*.

So the dominant welfare loss from centralisation is **Akerlof-style unravelling among the
people who never test at all** — a deadweight loss concentrated in transactions that never
happen, which is precisely why it is easy to miss: it is invisible in the data of people who
did test. It also reframes what this tool is. GenomeLens is not mainly "the private option
for people who would have tested anyway" — modelled this way it is an **access
intervention**, and in welfare terms access effects typically dominate quality effects.

*Why this matters methodologically:* the model converts a values argument ("privacy is
good") into **a falsifiable empirical question** — *is the annual breach hazard above
0.79%?* — and then the data answers it. Every parameter here is an assumption you can
change and re-run; none of it is a slogan.

**20 · Markov cohort model and budget impact — the two standard HEOR deliverables**

A static decision model cannot answer "what happens over 40 years?" or "can a payer afford
this?" Both canonical structures are implemented in [`markov_model.py`](markov_model.py).

**Markov state-transition cohort model.** A closed cohort moves between health states
(Well → Disease → Dead) on annual cycles under a validated transition matrix:

```math
\mathbf{s}_{t+1}=\mathbf{s}_t\mathbf{P},\qquad
p = 1-e^{-r\Delta t},\qquad
\text{ICER}=\frac{\Delta \text{Cost}}{\Delta \text{QALY}}
```

*In plain English:* track a group of people year by year. Each year some stay well, some
develop the disease, some die — of the disease or of anything else. Add up the costs and
quality-adjusted life-years each arm accrues, and the difference between the
genotype-guided arm and standard care gives the **ICER**, the number an HTA body reads.

Conventions a reviewer checks for, all implemented: **half-cycle correction** (states are
entered continuously, not on cycle boundaries), **rate→probability conversion**
`p = 1 − e^{−rΔt}` (never `r·Δt`, which can exceed 1), **age-dependent competing mortality**,
discounting of both costs and QALYs, and structural validation that the cohort is conserved
and death is absorbing.

**Budget impact analysis (BIA).** CEA asks *is it worth it?*; BIA asks *can we afford it?* —
and the conventions deliberately differ (ISPOR Task Force):

```math
\text{PMPM}_y=\frac{\text{Cost}_y^{\text{test}}+\text{Cost}_y^{\text{intervention}}-\text{Offsets}_y}{N_{\text{members}}\times 12}
```

*In plain English:* a payer doesn't care about 40-year QALYs — they care what hits next
year's budget. So BIA is short-horizon (1–5 years), scaled to a real plan population,
**undiscounted** (these are actual cash outlays), phases in **uptake** rather than assuming
instant adoption, and reports **per-member-per-month** — the metric that actually decides
formulary placement.

**Thresholds & sources.** `λ` = \$50k / \$100k / \$150k per healthy year (Neumann et al.,
*NEJM* 2014); `r` = 3% (Second Panel on Cost-Effectiveness in Health and Medicine);
cost-of-illness from the ADA, AHA, and Alzheimer's Association; drug–gene values from
published cost-effectiveness studies (Schackman 2008; Kazi 2014; Deenen 2016; CPIC);
health capital from Grossman (1972); option timing from Dixit & Pindyck (1994); EVPI from
Raiffa & Schlaifer (1961) and Claxton (1999); risk preferences from Arrow (1963) and Pratt
(1964); ascertainment correction from Begg (2002) and Gabai-Kapara (2014); winner's-curse
shrinkage from Zhong & Prentice (2008). Every figure is illustrative — a transparent model
you can inspect, not a clinical or financial guarantee.
