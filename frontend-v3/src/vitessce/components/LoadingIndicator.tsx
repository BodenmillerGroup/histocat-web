import { Spinner } from '@blueprintjs/core';
import React from 'react';

export default function LoadingIndicator() {
  return (
    <div className="loading-indicator-backdrop">
      <div className="loading-indicator-container">
        <Spinner />
      </div>
    </div>
  );
}
