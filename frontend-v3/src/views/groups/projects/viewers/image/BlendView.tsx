import styles from "./BlendView.module.scss";
import shallow from "zustand/shallow";
import React, { useState } from "react";
import { Intent } from "@blueprintjs/core";
import { useProjectsStore } from "modules/projects";
import { SuggestionView } from "components/SuggestionView";
import "vitessce/css/index.scss";
import Scatterplot from "vitessce/components/scatterplot/Scatterplot";
import { Status } from "vitessce/components/status";

export function BlendView() {
  const { activeAcquisition, selectedMetals, setSelectedMetals, getChannelStackImage } = useProjectsStore(
    (state) => ({
      activeAcquisition: state.getActiveAcquisition(),
      selectedMetals: state.selectedMetals,
      setSelectedMetals: state.setSelectedMetals,
      getChannelStackImage: state.getChannelStackImage,
    }),
    shallow
  );

  if (!activeAcquisition) {
    return <SuggestionView title="" content="Please select acquisition" intent={Intent.NONE} />;
  }

  const view = { target: [0, 0, 0], zoom: 0.75 };
  const mapping = "PCA";
  const cells = {
    1: { mappings: { [mapping]: [0, 0] } },
    2: { mappings: { [mapping]: [1, 1] } },
    3: { mappings: { [mapping]: [1, 2] } },
  };
  const selectedCellIds: any = [];
  const dimensions = { width: "800px", height: "400px", margin: "10px" };

  return (
    <div className={styles.container}>
      <span className={styles.toolbar}>
        <div className="vitessce-container vitessce-theme-light">
          <div className="card card-body bg-secondary" style={dimensions}>
            <Status info="Hello world" removeGridComponent={() => {}} />
          </div>
          <div className="card card-body bg-secondary" style={dimensions}>
            <Scatterplot
              uuid="my-vitessce-scatterplot"
              viewState={view}
              mapping={mapping}
              cells={cells}
              cellSelection={selectedCellIds}
              cellColors={null}
              setViewState={() => {}}
              updateStatus={(message: any) => {}}
              updateCellsSelection={(selectedIds: any) => {}}
              updateCellsHover={(hoverInfo: any) => {}}
              updateViewInfo={(viewInfo: any) => {}}
              clearPleaseWait={(layerName: any) => {}}
            />
          </div>
        </div>
      </span>
    </div>
  );
}
