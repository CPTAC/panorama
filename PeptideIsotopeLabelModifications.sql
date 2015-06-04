SELECT 
peptideisotopemodification.PeptideId.PeptideModifiedSequence,
peptideisotopemodification.IsotopeModId.Name ModificationName,
peptideisotopemodification.IsotopeModId.AminoAcid,
peptideisotopemodification.IsotopeModId.Terminus,
peptideisotopemodification.IsotopeModId.Formula,
peptideisotopemodification.IsotopeModId.Label13C,
peptideisotopemodification.IsotopeModId.Label15N,
peptideisotopemodification.IsotopeModId.Label18O,
peptideisotopemodification.IsotopeModId.Label2H,
runisotopemodification.IsotopeLabelId.Name IsotopeLabel,
(
CASE WHEN (peptideisotopemodification.IsotopeModId.Label13C = TRUE) 
THEN ( CASE WHEN(peptideisotopemodification.IsotopeModId.Label15N = TRUE OR 
                 peptideisotopemodification.IsotopeModId.Label18O = TRUE OR 
                 peptideisotopemodification.IsotopeModId.Label2H = TRUE) THEN('13C and ') ELSE('13C') END) 
ELSE ('') END || 
CASE WHEN (peptideisotopemodification.IsotopeModId.Label15N = TRUE) 
THEN ( CASE WHEN(peptideisotopemodification.IsotopeModId.Label18O = TRUE OR 
                 peptideisotopemodification.IsotopeModId.Label2H = TRUE) THEN('15N and ') ELSE('15N') END) 
ELSE ('') END || 
CASE WHEN (peptideisotopemodification.IsotopeModId.Label18O = TRUE) 
THEN ( CASE WHEN(peptideisotopemodification.IsotopeModId.Label2H = TRUE) THEN('18O and ') ELSE('18O') END)
ELSE ('') END ||
CASE WHEN (peptideisotopemodification.IsotopeModId.Label2H = TRUE) THEN ('2H') ELSE ('') END ||

CASE WHEN (peptideisotopemodification.IsotopeModId.Terminus = 'C') THEN (' at C-terminus') ELSE ('') END ||
CASE WHEN (peptideisotopemodification.IsotopeModId.Terminus = 'N') THEN (' at N-terminus') ELSE ('') END ||
CASE WHEN (peptideisotopemodification.IsotopeModId.AminoAcid IS NOT NULL) THEN (' ' || peptideisotopemodification.IsotopeModId.AminoAcid) ELSE ('') END
) AS IsotopeModification

FROM peptideisotopemodification, runisotopemodification
WHERE peptideisotopemodification.IsotopeModId = runisotopemodification.IsotopeModId