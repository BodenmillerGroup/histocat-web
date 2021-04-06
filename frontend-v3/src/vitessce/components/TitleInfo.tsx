import React from "react";
import classNames from "classnames";
import LoadingIndicator from "./LoadingIndicator";
import { Button } from "@blueprintjs/core";
import styles from "./TitleInfo.module.scss";

type TitleInfoProps = {
  title: string;
  info?: string;
  children?: any;
  isScroll?: boolean;
  isSpatial?: boolean;
  removeGridComponent(event: any): void;
  urls?: any[];
  isReady: boolean;
  options?: any;
  theme: string;
};

export default function TitleInfo(props: TitleInfoProps) {
  const { title, info, children, isScroll, isSpatial, removeGridComponent, isReady } = props;
  const childClassName = isScroll ? styles.scrollCard : isSpatial ? styles.spatialCard : styles.card;
  return (
    // d-flex without wrapping div is not always full height; I don't understand the root cause.
    <>
      <div className={classNames("title", styles.title)}>
        <span className={styles.titleContent}>
          {title}
        </span>
        <span className={styles.details}>
          <span className={styles.detailsContent}>
            {info}
            <Button onClick={removeGridComponent} icon="small-cross" minimal={true} />
          </span>
        </span>
      </div>
      <div className={childClassName}>
        {!isReady && <LoadingIndicator />}
        {children}
      </div>
    </>
    // "pl-2" only matters when the window is very narrow.
  );
}
