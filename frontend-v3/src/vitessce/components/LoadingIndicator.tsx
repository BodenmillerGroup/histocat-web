import { Spinner } from '@blueprintjs/core';
import React from 'react';
import styles from "./LoadingIndicator.module.scss";

export default function LoadingIndicator() {
  return (
    <div className={styles.backdrop}>
      <div className={styles.container}>
        <Spinner />
      </div>
    </div>
  );
}
