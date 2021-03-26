import styles from "./BlendView.module.scss";
import shallow from "zustand/shallow";
import React, { useMemo, useState } from "react";
import { Intent } from "@blueprintjs/core";
import { useProjectsStore } from "modules/projects";
import { SuggestionView } from "components/SuggestionView";
import "vitessce/css/index.scss";
import Scatterplot from "vitessce/components/scatterplot/Scatterplot";
import { Status } from "vitessce/components/status";
import { useSetViewConfig } from "vitessce/app/state/hooks";
import VitessceGrid from "vitessce/app/VitessceGrid";
import { getComponent } from "vitessce/app/component-registry";

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
    name: 'Linnarsson',
    version: '1.0.0',
    description: "SSSSS",
    public: true,
    datasets: [
      {
        uid: 'linnarsson-2018',
        name: 'Linnarsson 2018',
        description: `Linnarsson: SSSS`,
        files: [],
      },
    ],
    initStrategy: 'auto',
    coordinationSpace: {
      embeddingZoom: {
        PCA: 0,
        TSNE: 0.75,
      },
      embeddingType: {
        PCA: 'PCA',
        TSNE: 't-SNE',
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
      { component: 'spatial',
        coordinationScopes: {
          spatialZoom: 'A',
          spatialTargetX: 'A',
          spatialTargetY: 'A',
        },
        x: 2, y: 0, w: 4, h: 4 },
      { component: 'scatterplot',
        coordinationScopes: {
          embeddingType: 'PCA',
          embeddingZoom: 'PCA',
        },
        x: 6, y: 0, w: 3, h: 2 },
      { component: 'scatterplot',
        coordinationScopes: {
          embeddingType: 'TSNE',
          embeddingZoom: 'TSNE',
        },
        x: 6, y: 2, w: 3, h: 2 },
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
          <div className="card card-body bg-secondary" style={dimensions}>
            <VitessceGrid
              config={config}
              getComponent={getComponent}
              rowHeight={200}
              height={800}
              theme="dark"
            />
          </div>
        </div>
      </span>
    </div>
  );
}
