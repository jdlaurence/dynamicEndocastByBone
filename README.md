# dynamicEndocastByBone

Quantify the contribution of bones to oral cavity volume change. 

See "Volumetric analysis" methods in our Journal of Experimental Biology paper: 

[Go to instructions](https://github.com/jdlaurence/dynamicEndocastByBone/blob/main/instructions.md)

This method uses XROMM data and expands upon the dynamicEndocast function (https://www.pnas.org/content/112/28/8690). This package digitally freezes each bone in an animation relative to a reference bone, allowing measurement of the impact of the frozen bone (volume change in unaltered behavior - volume change when bone of interest is frozen). 

Freeze increments should be selected based upon frequency of behavior (we recommend roughly 10% of behavior duration). Researchers should additionally examine the "delta volume" curves at various freeze increments to check for oversmoothing.
