import styles from "./ProjectView.module.scss";
import React from "react";
import { TITLE_MAP, useLayoutsStore } from "modules/layouts";
import shallow from "zustand/shallow";
import { ViewId } from "modules/layouts/models";
import { getComponent } from "vitessce/app/component-registry";
import VitessceGrid from "vitessce/app/VitessceGrid";
import { useSetViewConfig } from "vitessce/app/state/hooks";

export function ProjectView() {
  const { setActiveNode, activeLayout } = useLayoutsStore(
    (state) => ({
      setActiveNode: state.setActiveNode,
      activeLayout: state.activeLayout,
    }),
    shallow
  );

  // const handleLayoutChange = (node: MosaicNode<ViewId> | null) => {
  //   setActiveNode(node);
  // };

  // const handleCreateNode = (): ViewId => {
  //   return "new";
  // };

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
        h: 12,
      },
      {
        component: "image",
        x: 2,
        y: 0,
        w: 8,
        h: 12,
      },
      {
        component: "channels",
        x: 10,
        y: 0,
        w: 2,
        h: 6,
      },
      {
        component: "settings",
        x: 10,
        y: 6,
        w: 2,
        h: 6,
      },
    ],
  };
  setViewConfig(config);

  return (
    <VitessceGrid config={config} getComponent={getComponent} height={1200} theme="dark" />
  );
}
