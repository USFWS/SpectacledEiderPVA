---
output: 
  pdf_document:
    latex_engine: xelatex
---
# Decision Memo — Reclassification Criterion for Uplisting the Alaska-Breeding Spectacled Eider Populations from Threatened to Endangered

**To:** Spectacled Eider Recovery Team; USFWS Northern Alaska Fish and Wildlife Field Office  

**From:** Erik E. Osnas, Biometrician, USFWS Migratory Bird Management  

**Re:** A coherent loss function for the loss-comparison criterion, and recommended decision framing  

**Status:** Draft for team discussion — not for distribution  

---

## 1. Purpose and the decision in front of us

Both Alaska breeding populations have declined, and the relevant classification question is now reclassification from threatened to **endangered**. The decline has been most dramatic on the YKD, but the ACP population also warrants evaluation against the criteria. The recovery plan's endangered criterion (USFWS 1996, p. 36) is an AND of two conditions:

> (1) the population is declining by $\ge$ 5%/year, judged by a loss comparison in which *the under-protection loss exceeds the over-protection loss*, calculated from trend data (5–15 yr) with loss functions where the loss when classifying is zero when r $\le$ -0.05 and the loss when not classifying is zero when r $\ge$ 0; **AND**
> (2) the minimum estimated population size is <3,000 breeding pairs for $\ge$ 1 year; **OR** <2,000 pairs in any one year (with a 2-year survey caveat).

**The abundance condition is met for both populations.** On the YKD, minimum estimated population size has been <3,000 pairs for more than one recent year. On the ACP, both the posterior mean and the lower 95% bound of estimated abundance in the last two years fall below the stricter <2,000-pair branch of the criterion: the posterior means are ≈3,747 and ≈3,801 individuals (≈1,873 and ≈1,901 pairs), and the lower 95% bounds are ≈2,532 and ≈2,434 individuals (≈1,266 and ≈1,217 pairs). The ACP therefore clears the <2,000-pair branch on both the central and the lower estimate. The abundance condition does not bear on either decision; in both cases the AND reduces to the **loss comparison** in condition (1). That comparison is the pivotal, and the contested, part of the decision — which is why getting it coherent is essential.

This memo (a) identifies the one quantity on which the loss comparison should rest, (b) shows why the recovery plan's original loss machinery is no longer a defensible way to compute it, (c) recommends how to close the decision in a way that is coherent, consistent with best available science, and transparent about the value judgments any classification necessarily embeds, and (d) addresses the specific objection that uplisting should not proceed because no management actions would follow from it.

---

## 2. The coherent loss quantity: the IPM 40-year quasi-extinction probability

The substance the recovery plan's loss function was always trying to capture is the **probability of reaching quasi-extinction (QE; <250 adults) within the management horizon.** The plan says so directly: "The cost of making decision errors (loss functions) is measured in terms of the probability of decreasing to under 250 adults in 50 years" (USFWS 1996, p. 30; Taylor et al. 1996). The stochastic random-walk model with drift (hereafter, drift model) in Appendix I/II was simply the tool available in 1996 for computing that probability.

That tool has been superseded. The 2026 Integrated Population Model–PVA computes the same quantity — P(QE within the horizon) — using the best available science: a density-dependent, age- and breeding-state-structured model that carries process variance, lead exposure, and winter sea-ice effects internally, and projects each draw forward from the current posterior state. The per-draw QE outcome in the IPM-PVA **is** the recovery plan's loss quantity, computed with the structure this population's dynamics require.

**I recommend the loss be defined as the IPM-PVA quasi-extinction probability, p = P(QE within 40 years).** The under-protection loss (failing to uplist a population genuinely at risk) attaches to this p directly, because p *is* the modeled risk the species faces if not given added protection.

For the two Alaska breeding populations, the IPM-PVA 40-year QE probabilities (Osnas 2026 draft, Tables 1 and 2) are:

| Scenario | YKD Pr(QE) by 2066 | ACP Pr(QE) by 2066 |
|---|---|---|
| Current Conditions | 0.05 | 0.11 |
| RCP 4.5 | 0.06–0.08 | 0.07 |
| RCP 8.5 | 0.09–0.11 | 0.11 |
| Continued lead (YKD only) | up to 0.16 | — |

(YKD scenarios vary by lead-deposition assumption as well as climate; the ACP model carries no lead term.) So the under-protection loss is **p ≈ 0.05 (YKD) and p ≈ 0.11 (ACP) under current conditions**, rising to 0.13–0.16 for the YKD under continued lead deposition combined with climate-driven ice loss. These are no longer trivial values. Note that the ACP current-conditions loss is the larger of the two, so the loss criterion is at least as strong for the ACP as for the YKD.

### What "loss" means here, and the one question the decision turns on

It is worth being explicit about what p is doing, because the word "loss" is carrying three distinct ideas that a decision should keep separate:

- **The outcome cost:** what it costs if the population actually reaches quasi-extinction — a large (not infinite) cost, since quasi-extinction at 250 birds is the point where only last-resort measures like captive breeding remain.
- **The probability of that outcome, p:** what the IPM computes — the chance the population reaches quasi-extinction within the horizon if it is not given added protection.
- **The decision loss, or regret:** what we actually minimize — not extinction itself, but the cost of getting the *classification* wrong. We regret under-protection if we fail to uplist a population that was genuinely headed for quasi-extinction; we regret over-protection if we uplist one that was not.

The recovery plan measures loss in units of quasi-extinction. Under that convention the under-protection loss is simply the expected number of quasi-extinctions, which for one population over one horizon equals p = P(QE). That is why p can stand in for "loss" — but it remains a probability, and it becomes a decision only when weighed against the cost of the opposite error.

That weighing is the one unavoidable value judgment, and it reduces to a single question a non-quantitative decision-maker can answer directly: **how many times worse is a quasi-extinction than an unnecessary regulatory burden?** That ratio is the cost ratio k used below. Everything else in the decision is science (the value of p) or arithmetic (the comparison); only k is a value judgment, and stating it openly — rather than burying it in a "symmetric" loss or an abundance threshold — is the central recommendation of this memo.

*Horizon note.* I use the 40-year QE horizon reported by the IPM-PVA, matching the horizon used in the 2021 SSA. I flag that the horizon has not been applied consistently across documents (the recovery plan and Taylor et al. (1996) reference 50 years; the IPM-PVA reports 40 and 74). I recommend the team adopt 40 years explicitly and state it as a decision convention, so the inconsistency is resolved going forward rather than carried silently.

---

## 3. Why the original loss machinery should not be used

Three features of the 1996 loss computation do not survive scrutiny. Two were already raised by the team in the Discussion of Dunham et al. (2021) and are properties of the recovery plan itself.

**(a) The drift-model surrogate omits the structure that governs this population.** In the 1996 plan, the drift model's process variance (σ) was a prior guess based on Common Eiders — two populations, one of which included an oil-spill die-off (Taylor et al. 1996). That particular deficiency has since been remedied: as of the 2021 SSA, σ is estimated directly from the spectacled eider data using the log-random-walk model, so a drift model fit today is no longer disqualified by a borrowed, cross-species σ. The deeper problem is structural. A simple drift model — however its σ is obtained — cannot represent the processes that actually drive this population's risk: density dependence, age and breeding-state structure, lead exposure, and the winter sea-ice effects on survival and breeding propensity. The IPM carries all of these internally and estimates the relevant variance structure from the data, while removing the need for the drift model's separate N₀ and growth-rate-grid choices. The case for the IPM is therefore not that the drift σ is untenable (it no longer is), but that the IPM estimates the quasi-extinction probability p directly, under the best available science for this population's dynamics, whereas the drift model is the wrong structure for computing the quantity the loss depends on.

**(b) The "zero loss when r $\ge$ 0 (or r $\le$ -0.05)" region convention double-handles risk.** The convention zeroes the loss based on where a point summary of r falls. That made sense only as a patch for a surrogate that could speak to nothing but r. Once the loss *is* the QE probability, the patch is not merely unnecessary — it is incoherent: a stable or slightly positive point estimate of r can still carry real, modeled quasi-extinction risk (from process variance, current low abundance, density dependence, and winter conditions), and the convention would zero out exactly that risk. The QE probability already integrates over r and everything else. Masking it by r counts the same information twice. **I recommend dropping the region convention and letting the QE probability stand as the whole loss.**

**(c) "Loss = 0 at r = 0" understates extinction risk.** As the team noted in Dunham et al. (2021), there is non-zero QE risk even when growth is non-negative, so a loss function pinned to zero at r = 0 does not reflect the actual risk that predicates the decision. Defining loss as p resolves this automatically.

---

## 4. Closing the decision: three coherent options

Defining the under-protection loss as p fixes the *loss*, but a classification decision is a *comparison*. The endangered rule compares under-protection (fail to uplist a population headed for QE) against over-protection (uplist one that is not). A coherent counterpart to p is needed for the over-protection side. Three options follow. All three define loss on the quantity we actually believe (the probability), and so all are more coherent than the recovery plan's reflected-r-curve.

**Option A — Explicit cost ratio (recommended).** Uplist if and only if (hereafter, iff) p·C_under > (1-p)·C_over, carrying a single interpretable ratio k = C_under/C_over for the relative cost of under- vs. over-protection. Equivalently, uplist iff p > p\* where the break-even threshold p\* = 1/(1+k). This is the cleanest expression of what the decision actually turns on, and it makes the one unavoidable value judgment — the relative cost of the two errors — explicit and inspectable.

**Option B — Threshold / Bayes rule (fallback).** Uplist iff p > p\* for an agreed acceptable QE probability p\*. This is the same rule as Option A with p\* stated directly rather than derived from k (p\* = 1/(1+k)); the two are interchangeable. It is the recommended fallback if the team prefers to discuss an acceptable-risk threshold rather than a cost ratio. Its known risk — that it hands decision-makers a dial they can set permissively — is mitigated by the break-even delivery in §5.

**Option C — Symmetric in probability (only if explicitly requested).** Define over-protection loss as the mirror of p about p\*, preserving the recovery plan's "equalize errors" intent but on p rather than r. I do not recommend this unless a decision-maker explicitly wants risk-neutrality: symmetry is not a neutral default, it is a specific 1:1 cost assumption (k = 1), and presenting it as the baseline would repeat the very move this memo is trying to correct. If a decision-maker affirmatively prefers equal weighting of the two errors, Option C is simply Option A with k = 1, and it can be adopted on that explicit basis.

### How these relate to the past methods, and their problems

The recovery plan and Taylor et al. (1996) implement a **risk-neutral, symmetric-around-r = 0** comparison — effectively Option C defined on r, with equalized errors. Cummings et al. (2018) classify this as the "putting species in the right bin" framing (their framing 1): compare a status assessment to thresholds, accounting for misclassification costs. My recommendation keeps that framing — appropriate and legally comfortable under the ESA's "best scientific data available" standard — while fixing two defects Cummings et al. identify as pervasive and largely unexamined in ESA practice:

- **The symmetric loss is itself a value choice, not a neutral default.** Equalizing over- and under-protection errors encodes a 1:1 cost ratio. It is no less a value judgment than any other ratio; it merely hides the judgment inside the word "symmetric" (Cummings et al. 2018; Dunham et al. 2021).
- **The abundance threshold embeds a hidden, and extreme, cost ratio.** The plan's "minimum estimated population size" uses the lower 95% bound of the estimate. As Dunham et al. (2021) show, a 2.5-percentile threshold implies, under a linear loss, that **overestimates of abundance are ~39× worse than underestimates** — an implied risk attitude far stronger than the symmetric-around-zero loss used elsewhere in the same criterion, and internally contradictory with it.

The combined force of these two points is the central argument for the team's deliberations with decision-makers: **refusing to state costs does not avoid a cost assumption; it only moves the assumption somewhere unexamined and unaccountable.** Every rule that converts p into an uplist/don't-uplist decision embeds a cost ratio whether or not anyone writes it down. The honest options make it visible.

---

## 5. Recommendation

**I recommend the explicit expected-cost framing (Option A), delivered operationally as a break-even / sensitivity analysis** so that a reluctant group is not forced to commit to an absolute cost. Concretely:

1. **Report p** — the IPM-PVA 40-year QE probability — as the under-protection loss: **0.05–0.16 for the YKD and 0.07–0.11 for the ACP, depending on scenario** (§2), carried across scenarios and reported under both ice-effect priors (§6) rather than as a single point.
2. **Report the break-even cost ratio** k\* at which the decision flips. Recall k is the answerable value question from §2: *how many times worse is a quasi-extinction than an unnecessary regulatory burden?* Uplisting is the lower-loss action whenever k exceeds (1-p)/p. Across the scenario range this break-even value runs, for the YKD, from **19** at the low end (p = 0.05) down to **≈5.3** at the high end (p = 0.16); for the ACP, from **≈13** (p = 0.07) to **≈8** (p = 0.11). In words: at the YKD low end, uplisting is favored unless the team judges a quasi-extinction to be *less than 19 times worse* than an unnecessary regulatory burden.
3. **Put the decision to the team and decision-makers as that single comparison:** *"How many times worse is a quasi-extinction of this population than the burden of protecting it unnecessarily? If the answer is more than [~5 to ~19, depending on scenario], uplisting is the lower-loss choice."* This replaces "pick a number" with a judgment people can actually make, and it makes the implied risk attitude required to decline uplisting explicit: declining to uplist means judging a quasi-extinction to be only a handful of times worse — or less — than an unnecessary regulatory burden.

**Fallback:** if the team prefers to discuss an acceptable-risk threshold rather than a cost ratio, use **Option B**, which is the identical rule with p\* stated directly. **Option C** should be used only if a decision-maker explicitly requests equal weighting (k = 1).

I recommend **against** reconstructing the drift-model loss, for the reasons in §3.

**A note on regulatory burden.** The regulatory burden of endangered status is the natural content of the over-protection cost against which k is weighed (§2) — an *error* cost, incurred when uplisting proves unwarranted. But regulatory burden is a socioeconomic cost, and the ESA directs that listing determinations be made "solely on the basis of the best scientific and commercial data available." It is therefore likely **inadmissible at the classification stage** — reserved for post-listing determinations (e.g., 4(d) rules, critical habitat) — the Cummings et al. (2018) framing-4 problem in concrete form.

This recommendation is consistent with the SSA's best-available-science posture: the recovery plan specified a simple drift model as a proxy for quasi-extinction risk; the IPM estimates that same risk directly and better, representing the density dependence, age and breeding structure, lead, and ice effects that the drift model omits; therefore the loss is evaluated using the IPM's quasi-extinction projection, and the decision rule is stated with its cost assumption explicit.

---

## 6. How "no available actions" relates to the decision

One consideration that may arise is that uplisting would promote no new management action, and so serves no purpose. It is worth locating where, if anywhere, this enters the decision structure.

The availability of management actions bears on the classification only through the cost of under-protection: "no actions would follow" is, in effect, the claim that failing to give the endangered classification costs nothing (C_under = 0). That is the only channel by which action-availability enters the loss comparison. Three points follow.

First, C_under = 0 is not generally true even where new recovery actions are limited. Endangered status changes Section 7 consultation standards (jeopardy determinations, incidental-take thresholds), can affect permitting and funding prioritization among listed species, strengthens the basis for the lead-shot and harvest measures the SSA identifies as pertinent, and shifts the baseline against which future actions are evaluated.

Second, under the break-even framing (§5), C_under = 0 makes the break-even cost ratio infinite: no quasi-extinction probability — not 0.16, not 0.5, not 1.0 — would ever justify uplisting. Declining to uplist on this basis is therefore equivalent to accepting certain quasi-extinction rather than changing the classification, a value statement made explicit once stated in these terms.

Third, the coherent home for a "no actions" concern is the timing of action, not the classification. Even where no action is available today, classification records the correct status, triggers reassessment, and positions the species for action when options emerge — Cummings et al.'s (2018) framing 2, "doing right by the species over time." The classification question ("what is the species' status?") and the action question ("what, if anything, do we do about it?") are distinct; the second does not answer the first.

---

## 7. Two issues the team should decide, and one caveat to state plainly

**Decision 1 — Scenario basis.** The choice of scenario is a judgment for the team, and I leave it open rather than prescribe one. If a soft lean is useful: RCP 4.5 is a reasonable central climate basis, as it is one of the two scenarios carried through the SSA and avoids anchoring on either the most optimistic or the most severe emissions path. The lead assumption is less clear to me, and I would defer to the team on which deposition scenario best represents current and expected conditions on the YKD. The choice matters: for the YKD the under-protection loss ranges from 0.05 to 0.16 across scenarios (break-even cost ratio ~19 down to ~5), and for the ACP from 0.07 to 0.11 (break-even ~13 down to ~8). Whatever basis is chosen, I recommend reporting the full range as a sensitivity envelope rather than a single point.

**Decision 2 — Horizon.** I recommend formally adopting the 40-year QE horizon and stating it as a convention (§2).

**Caveat — the winter-ice effect is now genuinely uncertain, and it affects the value of p.** The 2026 IPM-PVA finds the quadratic sea-ice effect on survival much weaker and less certain than in 2021 (posterior mean roughly one-third of the 2021 value; a substantial share of posterior draws are not concave-down in ice). Because p is the proportion of posterior trajectories reaching quasi-extinction, the modeling choice about how to represent the ice effect changes the value of p itself. In particular, the draft PVA notes that an informed prior forbidding survival-increasing-with-ice-deviation components would raise the estimated QE probability — a higher p, and therefore a lower break-even cost ratio, strengthening the case for uplisting. This is a point in favor of the QE-based loss, not against it: it surfaces the key scientific uncertainty at exactly the place it bears on the decision, rather than concealing it inside a drift σ. **I recommend reporting p under the current model alongside p recomputed under the alternative ice-effect prior, and treating the difference as a model-specification sensitivity.**

---

## 8. One-line summary for the SSA / 5-year review

> The endangered abundance condition is met for both Alaska breeding populations, so in each case the decision rests on the loss comparison. The coherent loss is the IPM-estimated 40-year quasi-extinction probability (0.05–0.16 for the YKD and 0.07–0.11 for the ACP, depending on scenario), which supersedes the recovery plan's drift-model loss and dispenses with its r-region convention. I recommend an explicit expected-cost comparison delivered as a break-even analysis, reducing the decision to one answerable question: *how many times worse is a quasi-extinction of this population than the burden of protecting it unnecessarily?* If more than roughly 5 to 19 times (depending on scenario), uplisting is the lower-loss choice. The objection that "no actions would follow" is an implicit claim that under-protection is costless, and conflates the classification decision with the separate question of what to do about it.

---

### Key sources relied on
- USFWS (1996), *Spectacled Eider Recovery Plan*, Recovery §A and Appendices I–II (loss functions, QE definition, abundance thresholds).
- Taylor, Wade, Stehn & Cochrane (1996), *Ecol. Appl.* 6:1077–1089 (original Bayesian loss/decision framework; symmetric, risk-neutral losses; σ from Common Eiders).
- Dunham, Osnas, Frost, Fischer & Grand (2021), *PLOS ONE* 16(7): e0253895. (critique of the recovery criteria: symmetric-loss and 2.5-percentile cost assumptions; redundancy of abundance and loss criteria; non-zero risk at r = 0).
- Cummings, Converse, Smith, Morey & Runge (2018), *Conserv. Biol.* 32:1246–1254 (implicit decision framing; misclassification-cost matrices; framing 1 vs. framings 2/4/5; the cost of leaving framing implicit).
- Osnas (2026 draft), *IPM-PVA for Spectacled Eiders: 2026 revisions* (40- and 74-year QE probabilities by scenario, Table 1; the revised, weaker, more uncertain ice effect).
- USFWS (2021), *SSA for the Spectacled Eider* (40-year QE horizon; Appendix A evaluation of the 1996 criteria).
