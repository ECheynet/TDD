# Mode shapes extraction by time domain decomposition (TDD)

[![View Mode shapes extraction by time domain decomposition (TDD) on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/52276-mode-shapes-extraction-by-time-domain-decomposition-tdd)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3817999.svg)](https://doi.org/10.5281/zenodo.3817999)


## Summary

The Time domain decomposition (TDD) [1] is an output-only method to extract mode shapes of a structure. Here, the modal damping ratios and modal displacements are in addition extracted using the functions presented in [6]. The TDD is similar to a more popular technique called Frequency-domain method (FDD) that was introduced by [2,3]. A good example of the FDD already exists on the Matlab File Exchange [4]. In a previous version, the present submission contained a function for the FDD. This function has been modified and moved to a new submission [5].


## Content

The submission contains:
- The function TDD.m: function to apply the TDD method.
- An example file Example1.m
- Acceleration data beamData.m (4 Mb)

## References

[1] Byeong Hwa Kim, Norris Stubbs, Taehyo Park, A new method to extract modal parameters using output-only responses, Journal of Sound and Vibration, Volume 282, Issues 1–2, 6 April 2005, Pages 215-230, ISSN 0022-460X, http://dx.doi.org/10.1016/j.jsv.2004.02.026.

[2] Brincker, R.; Zhang, L.; Andersen, P. (2001). "Modal identification of output-only systems using frequency domain decomposition". Smart Materials and Structures 10 (3): 441. doi:10.1088/0964-1726/10/3/303.

[3] BRINCKER, Rune, ZHANG, Lingmi, et ANDERSEN, P. Modal identification from ambient responses using frequency domain decomposition. In: Proc. of the 18*‘International Modal Analysis Conference (IMAC), San Antonio, Texas. 2000

[4] http://www.mathworks.com/matlabcentral/fileexchange/50988-frequency-domain-decomposition--fdd-

[5] https://se.mathworks.com/matlabcentral/fileexchange/57153-automated-frequency-domain-decomposition--afdd-

[6] https://se.mathworks.com/matlabcentral/fileexchange/55557-modal-parameters-identification-from-ambient-vibrations--sdof
