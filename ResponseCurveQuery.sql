SELECT
tci.TransitionId.PrecursorId.PeptideId.PeptideGroupId.Label AS Protein,
tci.TransitionId.PrecursorId.PeptideId.PeptideModifiedSequence AS PeptideModifiedSequence,
tci.TransitionId.PrecursorId.IsotopeLabelId.Name AS IsotopeLabel,
tci.TransitionId.PrecursorId.Charge AS PrecursorCharge,
tci.TransitionId.Charge AS ProductCharge,
CONCAT(tci.TransitionId.FragmentType, tci.TransitionId.FragmentOrdinal) AS FragmentIon,
tci.SampleFileId.ReplicateId.Name AS Replicate,
annot1.Value AS Concentration,
annot2.Value AS SampleGroup,
tci.Area,
tci.Background,
annot3.Value AS DoNotUse,
annot4.Value AS ISSpike,
--(SELECT annot1.Value FROM replicateannotation annot1 WHERE tci.SampleFileId.ReplicateId = annot1.replicateid AND annot1.Name='Concentration') AS Concentration,
--(SELECT annot2.Value FROM replicateannotation annot2 WHERE tci.SampleFileId.ReplicateId = annot2.replicateid AND annot2.Name='SampleGroup') AS SampleGroup,
--(SELECT annot3.Value FROM precursorchrominfoannotation annot3 WHERE tci.precursorChromInfoId = annot3.precursorChromInfoId AND annot3.Name='do not use') AS DoNotUse,

pepannot1.Value AS PeptideConcentrationIS,
pepannot2.Value AS PeptideConcentration,
annot5.Value AS MultiplicationFactor

FROM transitionchrominfo tci

LEFT JOIN peptideannotation pepannot1 ON (tci.TransitionId.PrecursorId.PeptideId = pepannot1.PeptideId AND pepannot1.Name='PeptideConcentrationIS')
LEFT JOIN peptideannotation pepannot2 ON (tci.TransitionId.PrecursorId.PeptideId = pepannot2.PeptideId AND pepannot2.Name='PeptideConcentration')
LEFT JOIN replicateannotation annot1 ON (tci.SampleFileId.ReplicateId = annot1.replicateid AND annot1.Name='Concentration')
LEFT JOIN replicateannotation annot2 ON (tci.SampleFileId.ReplicateId = annot2.replicateid AND annot2.Name='SampleGroup')
LEFT JOIN replicateannotation annot5 ON (tci.SampleFileId.ReplicateId = annot5.replicateid AND annot5.Name='MultiplicationFactor')
LEFT JOIN precursorchrominfoannotation annot3 ON (tci.precursorChromInfoId = annot3.precursorChromInfoId AND annot3.Name='do not use')
LEFT JOIN replicateannotation annot4 ON (tci.SampleFileId.ReplicateId = annot4.replicateid AND annot4.Name='IS Spike')