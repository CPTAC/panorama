SELECT
tci.TransitionId.PrecursorId.PeptideId.PeptideGroupId.Label AS Protein,
tci.TransitionId.PrecursorId.PeptideId.PeptideModifiedSequence AS PeptideModifiedSequence,
tci.TransitionId.PrecursorId.IsotopeLabelId.Name AS IsotopeLabel,
tci.TransitionId.PrecursorId.Charge AS PrecursorCharge,
tci.TransitionId.Charge AS ProductCharge,
CASE WHEN tci.TransitionId.NeutralLossMass > 0
THEN
(tci.TransitionId.FragmentType || tci.TransitionId.FragmentOrdinal || '-' || round(tci.TransitionId.NeutralLossMass,0))
ELSE
(tci.TransitionId.FragmentType || tci.TransitionId.FragmentOrdinal)
END
AS FragmentIon,
tci.SampleFileId.ReplicateId.Name AS Replicate,
annot1.Value AS Concentration,
annot2.Value AS SampleGroup,
annot3.Value AS ISSpike,
tci.Area,
tci.Background,
(SELECT annot3.Value
 FROM precursorchrominfoannotation annot3
WHERE annot3.PrecursorChrominfoId = tci.precursorChromInfoId
AND annot3.Name = 'do not use') AS DoNotUse,
pepannot1.Value AS PeptideConcentrationIS,
pepannot2.Value AS PeptideConcentration,
annot4.Value AS MultiplicationFactor
FROM transitionchrominfo tci
LEFT JOIN replicateannotation annot1 ON (tci.SampleFileId.ReplicateId = annot1.replicateid AND annot1.Name='Concentration')
LEFT JOIN replicateannotation annot2 ON (tci.SampleFileId.ReplicateId = annot2.replicateid AND annot2.Name='SampleGroup')
LEFT JOIN replicateannotation annot3 ON (tci.SampleFileId.ReplicateId = annot3.replicateid AND annot3.Name='IS Spike')
LEFT JOIN replicateannotation annot4 ON (tci.SampleFileId.ReplicateId = annot4.replicateid AND annot4.Name='MultiplicationFactor')
LEFT JOIN peptideannotation pepannot1 ON (tci.TransitionId.PrecursorId.PeptideId = pepannot1.PeptideId AND pepannot1.Name='PeptideConcentrationIS')
LEFT JOIN peptideannotation pepannot2 ON (tci.TransitionId.PrecursorId.PeptideId = pepannot2.PeptideId AND pepannot2.Name='PeptideConcentration')