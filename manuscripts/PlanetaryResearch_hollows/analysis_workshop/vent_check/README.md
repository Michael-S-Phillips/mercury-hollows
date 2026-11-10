# Vent position check (2 Nov 2026)

Catalog: data/data/MER_Vent_Locations_Classifications.csv (309 vents). Its "Lat" and "Lon" headers are swapped:
true latitude = column "Lon", true east longitude = column "Lat". Appears to be Pegg et al. (2021), ORDO dataset 12033876;
the compound fraction among good-image vents (77/109 = 71%) matches Pegg's thesis. Provenance: the file's macOS quarantine
tag records a Safari download on 1 Mar 2024 18:16 UTC; last modified 5 Mar 2024. Byte-identical copy in
Hollows_Research/python/data/. So the swap is in the file as downloaded, unless it was edited between 1 and 5 Mar 2024.

The PI's 2024 notebook data/von_mises_fit.ipynb reads the columns as labeled: cell 11 fits a von Mises 'vent longitude'
distribution to what are vent latitudes (folded into 0-62 and 284-360 deg), and cell 19 bins vents by the as-labeled
columns for the correlation matrices. Neither result is used in paper v5 or the SSS proposal.

| Check | Swapped reading | As labeled / null |
|---|---|---|
| IAU-named faculae (29) with a vent within 25 km (faculae_crossmatch.py) | 24 (median 9.5 km) | 0 within 50 km |
| Basemap contrast, mean DN within 8 km / 25-45 km (vent_imagery_check.py) | median 1.09, 72% > 1 | 1.00 (100 km shift); 1.01 as labeled |
| Catalog "In Crater" flag vs position inside a catalog crater (vent_assignment_check.py) | 85% agree (yes: 90% inside) | 50% agree after 100 km shift |
| Contact sheet (vent_offset_and_sheet.py, vent_sheet.png) | most vents on or beside irregular pits | — |

IAU centers checked against the USGS Gazetteer for Agwo, Nathair, Suge, Madu, Ngu, Bibilava, Coatl (all match the
Wikipedia list used). Five faculae (Abeeso, Ejo, Ibab, Nzoka, Thueban) sit at 0.0 km from a vent, so some IAU centers were
probably taken from the same mapping; the cross-match is not fully independent of the catalog.

Faculae with no vent within 50 km (vent_outlier_faculae.png): Suge (bright spot in SE Rachmaninoff; not in catalog),
Madu (no clear facula at this stretch; named 2023), Coatl (centered on a fresh ~10 km crater; may be impact-related),
Bibilava Faculae (100 km group; catalog vents to its SE are the Ejo/Ibab cluster). Omissions, not misplacements.

Test 8 paper specification on vent subsets (vent_assignment_check.json), odds ratio (95% CI):
all 4.88 (2.89-8.25); In Crater = yes 4.29 (2.49-7.39); not poor image 5.20 (2.93-9.23); Facula = yes 5.47 (3.14-9.55);
good image only 5.82 (2.70-12.55); In Crater = yes and not poor 4.45 (2.45-8.08); assignment within 0.8 R 3.90 (2.25-6.76).
