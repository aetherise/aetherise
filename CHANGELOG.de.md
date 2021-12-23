[English 🇺🇸️](CHANGELOG.md) / **Deutsch 🇩🇪️**

### Version 1.2.0
Hinzugefügt:
- Kreuzvalidierung `-loocv`
- Höhenlage einstellbar mit `-altitude`
- Datenblätter aus Cleveland 1927–1929
- Anpassungen für die Verarbeitung der Datenblätter aus Cleveland 1927–1929

Geändert:
- Verzeichnis Shankland1955 entfernt


### Version 1.1.0
Hinzugefügt:
- Filter `-abs_drift`
- Handbuch in Englisch

Geändert:
- Statistiken und Filter `-drift` beziehen sich nun auf einfach auf die mittlere Drift, nicht mehr auf die sogenannte berichtigte mittlere Drift.
- Fehlermeldungen
- Hilfetexte
- Handbuch für aetherise
- Bei der Aggregation der Spektren wird die Frequenz 0 nicht mehr ausgegeben.

Berichtigt:
- In `data_sheets.ods` und damit auch in `data_sheets.pdf` in den Spalten _sign_ und _sign correct_ falsche Zeichen entfernt.


### Version 1.0.0
Hinzugefügt:
- Datenreduzierung mittels DFT

Geändert:
- Die Zusatzstatistik (`-stats`) bei der Minimierung (`-aggregate fit`) zeigt nun die Normalverteilung der Residuen, statt der Chi-Quadrat-Werte.

Berichtigt:
- Fehlerbalken für Millers Algorithmus bei Verwendung des Schalters `-single`.
- Einige Datenblätter.
