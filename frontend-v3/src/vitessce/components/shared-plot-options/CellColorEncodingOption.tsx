import React from 'react';
import Select from '@material-ui/core/Select';
import TableCell from '@material-ui/core/TableCell';
import TableRow from '@material-ui/core/TableRow';
import { capitalize } from '../../utils';
import { useStyles } from './styles';

type CellColorEncodingOptionProps = {
  observationsLabel: string;
  cellColorEncoding: any;
  setCellColorEncoding: any;
}

export default function CellColorEncodingOption(props: CellColorEncodingOptionProps) {
  const {
    observationsLabel,
    cellColorEncoding,
    setCellColorEncoding,
  } = props;

  const classes = useStyles();

  const observationsLabelNice = capitalize(observationsLabel);

  function handleColorEncodingChange(event: any) {
    setCellColorEncoding(event.target.value);
  }

  return (
    <TableRow>
      {/*<TableCell className={classes.labelCell} htmlFor="cell-color-encoding-select">*/}
      <TableCell className={classes.labelCell}>
        {observationsLabelNice} Color Encoding
      </TableCell>
      <TableCell className={classes.inputCell}>
        <Select
          native
          className={classes.select}
          value={cellColorEncoding}
          onChange={handleColorEncodingChange}
          inputProps={{
            id: 'cell-color-encoding-select',
          }}
        >
          <option value="cellSetSelection">Cell Sets</option>
          <option value="geneSelection">Gene Expression</option>
        </Select>
      </TableCell>
    </TableRow>
  );
}
