# Case Study: Using an AI Assistant to Debug a Complex Statistical Model

## Purpose

This document summarizes a working session in which an analyst used an AI assistant (Claude) to debug a complex ecological statistical model. It is intended both for leadership evaluating AI use cases and for technical staff who may recognize the specific methods. The goal is to show *how* the collaboration worked — what the AI contributed, what the analyst contributed, and why the combination was effective — rather than to teach the underlying statistics.

## Background (in plain terms)

The analyst was building a model to estimate the survival and detection of a bird population from years of banding and resighting records. The technical form is a **Bayesian multi-state Cormack-Jolly-Seber (CJS) mark-recapture model**, fit using JAGS software, with the added complication that some years had **no marking or resighting effort** ("gap years"). The model tracks animals through five life stages (duckling through breeding adult), three of which can never be directly observed. The end goal is to fold this survival model into a larger **integrated population model (IPM)**, so both correctness and computational speed mattered.

A key design decision shaped the whole effort. The data can be represented two ways. The intuitive way tracks each individual animal's full detection history one at a time; it is easier to write and reason about, but slow to fit — in this case, on the order of hours. The alternative, the **m-array format**, collapses those individual histories into summary counts of when animals released in one year were next seen, and fits in minutes. That speed advantage was essential, because the model will ultimately run inside a much larger integrated model where the cost compounds. The trade-off is that the m-array is considerably more complex to code correctly and less forgiving of small errors — several of the bugs below stemmed directly from that complexity. The analyst deliberately accepted the harder-to-write format to gain the computational speed the larger model requires.

Crucially, models like this are validated by *simulation*: the analyst generates fake data from known "true" values, fits the model, and checks whether it recovers those values. When it doesn't, something is wrong — and the error could be in the model code, the simulation, or the way data are passed between them. Isolating which is the hard part.

## How the collaboration was structured

The defining feature of this session was a disciplined, layered strategy: validate the simplest version of the model first, then add one layer of complexity at a time, confirming correctness at each step before proceeding. This approach was reinforced throughout by both parties and was decisive, because the model contained *multiple independent bugs that masked one another*. Testing everything at once would have made them nearly impossible to separate.

The analyst drove the investigation. They ran every model fit, supplied the actual code, reported precise symptoms, and — critically — **repeatedly corrected the AI when its hypotheses were wrong**. The AI contributed structural knowledge of the method, hypotheses about likely causes, hand-calculations to check specific numbers, and diagnostic tests to distinguish competing explanations.

## The bugs found, and how

Over the session, the team identified and fixed five distinct problems:

1. **A data-length error.** An indicator vector meant to flag effort years was the wrong length, so gap years were never actually being switched off. *The analyst found this independently* while inspecting the data being passed to the model.

2. **A frozen time index in the model's core calculation.** A survival-accumulation loop reused the wrong year's parameters. The AI identified this by reading the code; it was subtle because it was numerically silent unless survival varied across years — which explained a small bias the analyst had noticed earlier.

3. **A "fall-through" logic error in the simulation.** A chain of conditional statements let simulated birds advance through multiple life stages in a single year, collecting extra detection and survival chances. This was the source of a persistent, puzzling bias in duckling survival. The team localized it by computing, by hand, the rate at which simulated ducklings *should* have been resighted (~0.032) versus the rate observed (~0.06) — a clean 2× discrepancy that pointed directly at the mechanism. *The analyst confirmed the fix* by re-running the check and recovering the expected value.

4. **Mishandled data when imposing gap years.** The analyst had zeroed out sightings in the finished data summary, which silently deleted birds instead of reclassifying them as "never seen again," biasing survival downward. *The analyst proposed the better approach* — imposing the gaps earlier in the pipeline — which also happened to mirror how real field data naturally arrive.

5. **An "off-by-one" alignment error.** The years switched off in the data did not line up with the years the model treated as gaps, because of a subtle difference in how release years and resighting years are indexed. The AI traced the exact mismatch; the fix aligned the two and resolved a final error.

After these fixes, the model correctly recovered all known values across every scenario tested, including the difficult cases involving gap years.

## A real phenomenon, not a bug

A final anomaly proved instructive. Years with very *low* (but nonzero) detection produced a large bias in the *following* year's detection estimate. The analyst suspected another bug. Investigation showed instead that this was a **genuine statistical property** of the data: when almost no animals are detected in a year, the information needed to separate survival from detection nearly disappears, and the two become confounded.

The team distinguished "real phenomenon" from "bug" through targeted tests. One was the analyst's suggestion to *reverse* the manipulation, making detection unusually **high** in some years, which produced the mirror-image bias pattern a static coding error could not plausibly cause. The decisive test was a simulation study comparing the fast m-array model against a slower, independently coded individual-history version of the same model: across replicates, both produced essentially identical patterns, showing the effect is shared between the two implementations and therefore not an artifact of the m-array structure. The precise statistical mechanism remains unclear. The two-model agreement rules out an error in the likelihood machinery — the code that differs between the m-array and individual-history implementations — so the team is confident there is no likelihood-code bug. It does not, however, eliminate a shared error in the parts both models have in common: the state-transition (`psi`) specification and the data simulation. Those possibilities have not yet been directly tested (a single-state CJS replication would do so). On present evidence a genuine statistical property is the most likely explanation, but the finding is best treated as a well-localized limitation to plan around in the larger model rather than a fully resolved one.

## What made this effective — and its limits

Several themes stand out for anyone considering similar use:

- **The human stayed in control and skeptical.** The AI was wrong more than once — at one point it confidently recommended a change that the analyst tested and found to break the model, and it retracted the advice. The analyst's willingness to run checks and push back was essential; treating the AI's output as authoritative would have introduced errors rather than removed them.

- **Concrete verification beat argument.** Progress was fastest when a specific number was computed and compared against a hand calculation or known truth, rather than reasoning about code in the abstract. The AI's ability to perform these checks on demand was valuable, but only because the analyst grounded them in real outputs.

- **Domain expertise was shared, not delegated.** The analyst supplied the biology, the actual code, and the judgment about what was plausible; the AI supplied method structure, pattern-recognition across a long debugging arc, and tireless bookkeeping of a dozen interacting indices. Neither party could have reached the result as efficiently alone.

- **The AI's main risks were overconfidence and plausible-but-wrong hypotheses.** These were caught only because the analyst verified independently. A less experienced user, or one inclined to trust the tool, could have been led astray.

## Bottom line

The session resolved five genuine bugs and correctly diagnosed one non-bug in a sophisticated statistical model, producing a validated result the analyst can build on. The value came not from the AI solving the problem autonomously, but from a tight loop between an expert analyst running tests and making judgments and an assistant contributing structure, hypotheses, and verification. The outcome depended on the analyst's expertise and skepticism as much as on the tool's capabilities.
