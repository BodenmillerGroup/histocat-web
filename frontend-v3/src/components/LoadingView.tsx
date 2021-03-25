import { NonIdealState, Spinner } from "@blueprintjs/core";
import styles from "./LoadingView.module.scss";

type LoadingViewProps = {
  title?: string;
  description?: string;
};

export function LoadingView(props: LoadingViewProps) {
  const { title = "Loading...", description = "Please wait while data is loading" } = props;
  return (
    <div className={styles.container}>
      <NonIdealState title={title} icon={<Spinner />} description={description} />
    </div>
  );
}
