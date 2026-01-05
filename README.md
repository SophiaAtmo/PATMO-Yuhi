# PATMO-Yuhi
[DOI Link](https://doi.org/10.5194/acp-2022-68)

v1 - From Seba

v2 - Adjust δ34S effective source signature in Manuscript 

v3 - Emission data consider NO MIF (δ33S=0.515*δ34S, δ36S=1.9*δ34S)

	COS	723 molec/cm3/s		δ34S=10.5 ‰	NO MIF
	
	CS2	394.1 molec/cm3/s		δ34S=10.4 ‰	NO MIF
	
	H2S	11011 molec/cm3/s	δ34S=1 ‰		NO MIF
	
	SO2	59081 molec/cm3/s	δ34S=5 ‰		NO MIF
	
	DMS	26406 molec/cm3/s	δ34S=20 ‰		NO MIF
	
        Dry deposition applies original manuscript data
		
	SO2	7.52025E-7 → 7.15E-7 s-1
	
	CS2 	29.97E-9 → 3.00E-9 s-1
	
        Simulation time  20 years → 60 years (COS reaches steady-state)
		
v4 - Compatible with Python 3 PATMO

        add input file
		
        add CS2 photo-oxidation
		
        update k
		
        shield reverse reaction
		
        COS emission parameter: 0.8 --> v4.1
