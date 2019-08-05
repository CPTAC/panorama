SELECT
tci.TransitionId.PrecursorId.PeptideId.PeptideGroupId.Label AS Protein,
tci.TransitionId.PrecursorId.PeptideId.PeptideModifiedSequence AS PeptideModifiedSequence,
tci.TransitionId.PrecursorId.IsotopeLabelId.Name AS IsotopeLabel,
tci.TransitionId.PrecursorId.Charge AS PrecursorCharge,
tci.TransitionId.Charge AS ProductCharge,

CASE WHEN tci.TransitionId.NeutralLossMass > 0
THEN
(tci.TransitionId.FragmentType || tci.TransitionId.FragmentOrdinal || '-' || ceiling(tci.TransitionId.NeutralLossMass))
ELSE
(tci.TransitionId.FragmentType || tci.TransitionId.FragmentOrdinal)
END
AS FragmentIon,


--CONCAT(tci.TransitionId.FragmentType, tci.TransitionId.FragmentOrdinal) AS FragmentIon,
tci.Area,
tci.SampleFileId.ReplicateId.Name AS ReplicateName,
(SELECT rannot2.Value FROM replicateannotation rannot2 WHERE tci.SampleFileId.replicateId = rannot2.replicateid AND rannot2.Name='Replicate Number') AS ReplicateNumber,
(SELECT rannot2.Value FROM replicateannotation rannot2 WHERE tci.SampleFileId.replicateId = rannot2.replicateid AND rannot2.Name='Temperature') AS Temperature,
--rannot2.Value AS Temperature
(SELECT rannot2.Value FROM replicateannotation rannot2 WHERE tci.SampleFileId.replicateId = rannot2.replicateid AND rannot2.Name='Time') AS `Time`,
--rannot2.Value AS Time
(SELECT rannot2.Value FROM replicateannotation rannot2 WHERE tci.SampleFileId.replicateId = rannot2.replicateid AND rannot2.Name='Time Units') AS TimeUnits,
--rannot2.Value AS TimeUnits
(SELECT rannot2.Value FROM replicateannotation rannot2 WHERE tci.SampleFileId.replicateId = rannot2.replicateid AND rannot2.Name='Freeze Thaw Cycles') AS FreezeThawCycles,
--rannot2.Value AS TimeUnits
(SELECT rannot3.Value FROM replicateannotation rannot3 WHERE tci.SampleFileId.replicateId = rannot3.replicateid AND rannot3.Name='Exp4 SampleGroup') AS Exp4SampleGroup,
--rannot3.Value AS SampleGroup,
--rannot4.Value AS DilutionFactor,
--pepannot.Value AS AnalyteConcentration

FROM transitionchrominfo tci

--INNER JOIN replicateannotation rannot1 ON (tci.SampleFileId.replicateId = rannot1.replicateid AND rannot1.Name='Replicate')
--INNER JOIN replicateannotation rannot2 ON (tci.SampleFileId.replicateId = rannot2.replicateid AND rannot2.Name='Concentration')
--INNER JOIN replicateannotation rannot3 ON (tci.SampleFileId.replicateId = rannot3.replicateid AND rannot3.Name='SampleGroup')
--LEFT JOIN replicateannotation rannot4 ON (tci.SampleFileId.replicateId = rannot4.replicateid AND rannot4.Name='DilutionFactor')
--LEFT JOIN peptideannotation pepannot ON (tci.TransitionId.PrecursorId.PeptideId = pepannot.peptideId AND pepannot.Name='Analyte Concentration')
