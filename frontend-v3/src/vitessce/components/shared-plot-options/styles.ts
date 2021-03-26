import { makeStyles } from "@material-ui/core/styles";

export const useStyles = makeStyles((theme) => ({
  box: {
    boxSizing: "border-box",
  },
  checkbox: {
    padding: "3px",
    color: (theme.palette as any).primaryForeground,
    "&:checked": {
      color: (theme.palette as any).primaryForeground,
    },
    "& input": {
      height: "100%",
    },
  },
  slider: {
    color: (theme.palette as any).primaryForeground,
    minWidth: "60px",
    padding: "10px 0 6px 0",
  },
  sliderValueLabel: {
    "& span": {
      "& span": {
        color: (theme.palette as any).primaryBackground,
      },
    },
  },
  tableContainer: {
    overflow: "hidden",
  },
  labelCell: {
    padding: "2px 8px 2px 16px",
  },
  inputCell: {
    padding: "2px 16px 2px 8px",
    overflow: "visible",
  },
  select: {
    "& select": {
      fontSize: ".875rem",
    },
  },
}));
