# Bayesian Causal Estimation with Stan.
## StanCon 2020 Presentation Repo

This folder contains slides (StanCon_BayesianCausal.pdf) and codes (code/) running synthetic example discussed in the talk.

Collapsed slides without LaTeX \pause is here /collapsed_slides/StanCon_BayesianCausal.pdf

This talk was based on the following paper:
https://arxiv.org/pdf/2004.07375.pdf

Feel free to cite this paper when sharing contents from this talk.

## Errata
Note original slides used during presentation is in old/

Aug. 16 2020:

	- When computing conditional causal effects, we want to integrate over P(L | V=v), not the marginal P(L). Slides and code are adjusted to reflect this. The original slides implicitly assumed the L and V are independent, therefore P(L) = P(L | V). It works fine because this is how the toy data was simulated to begin with, but is not general enough.
	- The sensitivity slide should be subtracting Delta off from the risk difference, NOT the expected risk under each treatment. This was corrected in the slide. 
	- The term pi wasn't defined in original slide, it is the conditional propensity score. I've added the definition to current slides.
  
