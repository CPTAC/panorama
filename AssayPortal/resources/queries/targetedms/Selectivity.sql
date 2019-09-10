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
tci.Area AS Area,
tci.SampleFileId.ReplicateId.Name AS ReplicateName,
(SELECT rannot1.Value FROM replicateannotation rannot1 WHERE tci.SampleFileId.replicateId = rannot1.replicateid AND rannot1.Name='Replicate Number') AS ReplicateNumber,
--rannot1.Value AS Replicate,

(SELECT rannot3.Value FROM replicateannotation rannot3 WHERE tci.SampleFileId.replicateId = rannot3.replicateid AND rannot3.Name='Exp3 SampleGroup') AS Exp3SampleGroup,
tci.SampleFileId.ReplicateId.SampleType AS SampleType,
--rannot3.Value AS SampleGroup,
--rannot4.Value AS DilutionFactor,
--pepannot.Value AS AnalyteConcentration
tci.SampleFileId.ReplicateId.AnalyteConcentration AS AnalyteConcentration,

FROM transitionchrominfo tci

--INNER JOIN replicateannotation rannot1 ON (tci.SampleFileId.replicateId = rannot1.replicateid AND rannot1.Name='Replicate')
--INNER JOIN replicateannotation rannot2 ON (tci.SampleFileId.replicateId = rannot2.replicateid AND rannot2.Name='Concentration')
--INNER JOIN replicateannotation rannot3 ON (tci.SampleFileId.replicateId = rannot3.replicateid AND rannot3.Name='SampleGroup')
--LEFT JOIN replicateannotation rannot4 ON (tci.SampleFileId.replicateId = rannot4.replicateid AND rannot4.Name='DilutionFactor')
--LEFT JOIN peptideannotation pepannot ON (tci.TransitionId.PrecursorId.PeptideId = pepannot.peptideId AND pepannot.Name='Analyte Concentration')
