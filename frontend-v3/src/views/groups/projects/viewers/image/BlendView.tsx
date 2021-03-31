import styles from "./BlendView.module.scss";
import shallow from "zustand/shallow";
import React, { useMemo, useState } from "react";
import { Intent } from "@blueprintjs/core";
import { useProjectsStore } from "modules/projects";
import { SuggestionView } from "components/SuggestionView";
import "vitessce/css/index.scss";
import { Status } from "vitessce/components/status";
import { useSetViewConfig } from "vitessce/app/state/hooks";

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
  const setViewConfig = useSetViewConfig();
  const config = {
    name: "Linnarsson",
    version: "1.0.0",
    description: "SSSSS",
    public: true,
    datasets: [
      {
        uid: "linnarsson-2018",
        name: "Linnarsson 2018",
        description: `Linnarsson: SSSS`,
        files: [],
      },
    ],
    initStrategy: "auto",
    coordinationSpace: {
      embeddingZoom: {
        PCA: 0,
      },
      embeddingType: {
        PCA: "PCA",
      },
      spatialZoom: {
        A: -5.5,
      },
      spatialTargetX: {
        A: 16000,
      },
      spatialTargetY: {
        A: 20000,
      },
    },
    layout: [
      {
        component: "slides",
        x: 0,
        y: 0,
        w: 2,
        h: 3,
      },
      {
        component: "channels",
        x: 2,
        y: 3,
        w: 2,
        h: 3,
      },
    ],
  };
  setViewConfig(config);

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
        </div>
      </span>
    </div>
  );
}
