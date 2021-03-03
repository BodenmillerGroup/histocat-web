import styles from "./SuggestionView.module.scss";
import { Callout, Intent } from "@blueprintjs/core";

type SuggestionViewProps = {
  title: string;
  content: string;
  intent: Intent;
};

export function SuggestionView(props: SuggestionViewProps) {
  return (
    <div className={styles.container}>
      <Callout title={props.title} intent={props.intent}>
        {props.content}
      </Callout>
    </div>
  );
}
