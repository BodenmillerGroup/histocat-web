import { NonIdealState, Spinner } from "@blueprintjs/core";
import styles from "./LoadingView.module.scss";

export function LoadingView() {
  return (
    <NonIdealState
      title="Loading..."
      icon={<Spinner />}
      description="Please wait while data is loading."
      className={styles.container}
    />
  );
}
