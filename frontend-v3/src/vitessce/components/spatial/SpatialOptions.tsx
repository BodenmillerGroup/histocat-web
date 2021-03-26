import React from 'react';
import OptionsContainer from '../shared-plot-options/OptionsContainer';
import CellColorEncodingOption from '../shared-plot-options/CellColorEncodingOption';

type SpatialOptionsProps = {
  observationsLabel: string;
  cellColorEncoding: any;
  setCellColorEncoding: any;
}

export default function SpatialOptions(props: SpatialOptionsProps) {
  const {
    observationsLabel,
    cellColorEncoding,
    setCellColorEncoding,
  } = props;

  return (
    <OptionsContainer>
      <CellColorEncodingOption
        observationsLabel={observationsLabel}
        cellColorEncoding={cellColorEncoding}
        setCellColorEncoding={setCellColorEncoding}
      />
    </OptionsContainer>
  );
}
