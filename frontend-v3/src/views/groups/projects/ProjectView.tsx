import styles from "./ProjectView.module.scss";
import React from "react";
import "react-mosaic-component/react-mosaic-component.css";
import { TITLE_MAP, useLayoutsStore } from "modules/layouts";
import shallow from "zustand/shallow";
import { ViewId } from "modules/layouts/models";
import { getComponent } from "vitessce/app/component-registry";
import VitessceGrid from "vitessce/app/VitessceGrid";
import { useSetViewConfig } from "vitessce/app/state/hooks";
import "vitessce/css/index.scss";

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
        h: 24,
      },
      {
        component: "channels",
        x: 10,
        y: 0,
        w: 6,
        h: 12,
      },
      {
        component: "settings",
        x: 10,
        y: 12,
        w: 2,
        h: 12,
      },
    ],
  };
  setViewConfig(config);

  return (
    <VitessceGrid config={config} getComponent={getComponent} height={1400} theme="dark" />
    // <div className={styles.container}>
    //   <Mosaic<ViewId>
    //     resize={{
    //       minimumPaneSizePercentage: 15,
    //     }}
    //     initialValue={activeLayout.node}
    //     onChange={handleLayoutChange}
    //     className={classNames("mosaic-blueprint-theme", Classes.DARK)}
    //     renderTile={(id, path) => (
    //       <MosaicWindow<ViewId> path={path} title={TITLE_MAP[id]}>
    //         {
    //           {
    //             slides: <SlidesView />,
    //             image: <BlendView />,
    //             channels: <ChannelsView />,
    //             settings: <ChannelsSettingsView />,
    //             new: <SettingsView />,
    //           }[id]
    //         }
    //       </MosaicWindow>
    //     )}
    //   />
    // </div>
  );
}
