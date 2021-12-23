**English 🇺🇸️** / [Deutsch 🇩🇪️](CHANGELOG.de.md)

### Version 1.2.0
Added:
- Cross-validation `-loocv`
- Option `-altitude`
- Data sheets from Cleveland 1927–1929
- Adjustments to process the data sheets from Cleveland 1927–1929

Changed:
- Removed directory Shankland1955


### Version 1.1.0
Added:
- Filter `-abs_drift`
- Manual in English

Changed:
- Statistics and filter `-drift` now simply refer to the mean drift, not to the so called corrected mean drift.
- Error messages
- Help text
- Manual of aetherise
- Removed the frequency 0 from the aggregated spectra

Fixed:
Berichtigt:
- In `data_sheets.ods` and thus also in `data_sheets.pdf` some wrong characters were removed in the columns _sign_ and _sign correct_.


### Version 1.0.0
Added:
- Data reduction using DFT

Changed:
- The additional statistics (`-stats`) at minimization (`-aggregate fit`) now shows the normal distribution of the residuals, instead of the chi-squared values.

Fixed:
- Error bars for Miller's algorithm if the option `-single` was used.
- Some data sheets.
