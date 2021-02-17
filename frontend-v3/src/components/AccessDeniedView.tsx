import styles from "./AccessDeniedView.module.scss";
import { Callout } from "@blueprintjs/core";

type AccessDeniedViewProps = {
  title: string;
  content: string;
};

export function AccessDeniedView(props: AccessDeniedViewProps) {
  return (
    <div className={styles.container}>
      <Callout title={props.title} intent="danger">
        {props.content}
      </Callout>
    </div>
  );
}
