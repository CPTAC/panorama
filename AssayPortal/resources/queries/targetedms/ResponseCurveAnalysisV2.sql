SELECT tci.TransitionId.PrecursorId.PeptideId.PeptideGroupId.Label          AS Protein,
       tci.TransitionId.PrecursorId.PeptideId.PeptideModifiedSequence       AS PeptideModifiedSequence,
       tci.TransitionId.PrecursorId.IsotopeLabelId.Name                     AS IsotopeLabel,
       tci.TransitionId.PrecursorId.Charge                                  AS PrecursorCharge,
       tci.TransitionId.Charge                                              AS ProductCharge,
       CASE
           WHEN tci.TransitionId.NeutralLossMass > 0
               THEN
               (tci.TransitionId.FragmentType || tci.TransitionId.FragmentOrdinal || '-' ||
                round(tci.TransitionId.NeutralLossMass, 0))
           ELSE
               (tci.TransitionId.FragmentType || tci.TransitionId.FragmentOrdinal)
           END
                                                                            AS FragmentIon,
       tci.Area,
       tci.Background,
       tci.SampleFileId.ReplicateId.Name                                    AS ReplicateName,
       (SELECT rannot1.Value
        FROM replicateannotation rannot1
        WHERE tci.SampleFileId.replicateId = rannot1.replicateid
          AND rannot1.Name = 'Replicate Number')                            AS ReplicateNumber,
       annot2.Value                                                         AS SampleGroup,
       tci.SampleFileId.ReplicateId.AnalyteConcentration                    AS AnalyteConcentration,
       tci.TransitionId.PrecursorId.PeptideId.InternalStandardConcentration AS InternalStandardConcentration,
       tci.TransitionId.PrecursorId.PeptideId.ConcentrationMultiplier       AS ConcentrationMultiplier,
       tci.SampleFileId.ReplicateId.SampleType                              AS SampleType,
       annot3.Value                                                         AS DoNotUse,

FROM transitionchrominfo tci
         LEFT JOIN replicateannotation annot2
                   ON (tci.SampleFileId.ReplicateId = annot2.replicateid AND annot2.Name = 'SampleGroup')
         LEFT JOIN precursorchrominfoannotation annot3
                   ON (tci.precursorChromInfoId = annot3.precursorChromInfoId AND annot3.Name = 'do not use')
