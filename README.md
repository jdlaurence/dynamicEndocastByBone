# dynamicEndocastByBone

_**Quantify the contribution of individual bones to oral cavity volume change in [XROMM](https://www.xromm.org/) animations.**_

**Publication:** [Whitlow et al. (2022)](https://journals.biologists.com/jeb/article-abstract/doi/10.1242/jeb.243283/273979/Suction-feeding-biomechanics-of-Polypterus-bichir?redirectedFrom=fulltext) -- See _Volumetric Analysis_ section of paper for a detailed explanation. 

**Description:** In short, `dynamicEndocastByBone` is a MATLAB function that performs a rolling freeze (for a user-specified duration) of each individual bone in an XROMM animation relative to a reference bone (e.g., neurocranium). The difference in endocast volume between the frozen and unfrozen animation is the impact of that bone's motion on endocast volume, at that time-point. The impact of freezing a given bone relative to the sum of all bones' impacts is that bone's **relative contribution to volume change (RCVC)**.

**Credits:** Concieved of by Katie Whitlow and J.D. Laurence-Chasen. Code written and maintained by J.D. Laurence-Chasen. Based on the original dynamicEndocast method by Ariel Camp, detailed in [Camp et al. (2015)](https://www.pnas.org/content/112/28/8690).

**[Go to instructions](https://github.com/jdlaurence/dynamicEndocastByBone/blob/main/instructions.md)**

![RCVC](https://user-images.githubusercontent.com/53494838/149544056-bbe0d0e4-7e69-44cc-bb6b-6d34200e7941.png)


