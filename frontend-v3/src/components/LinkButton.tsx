import { Button } from "@blueprintjs/core";
import { Link } from "react-router-dom";

type LinkButtonProps = {
  to: string;
  text: string;
  minimal?: boolean;
  icon?: any;
};

export function LinkButton(props: LinkButtonProps) {
  return (
    <Link to={props.to} style={{ textDecoration: "none" }}>
      <Button text={props.text} minimal={props.minimal} icon={props.icon} />
    </Link>
  );
}
