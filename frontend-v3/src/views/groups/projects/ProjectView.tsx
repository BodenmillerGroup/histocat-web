import styles from "./ProjectView.module.scss";
import React from "react";
import "react-mosaic-component/react-mosaic-component.css";
import { Mosaic, MosaicNode, MosaicWindow } from "react-mosaic-component";
import classNames from "classnames";
import { Classes } from "@blueprintjs/core";
import { TITLE_MAP, useLayoutsStore } from "modules/layouts";
import shallow from "zustand/shallow";
import { ViewId } from "modules/layouts/models";
import { BlendView } from "./viewers/image/BlendView";
import { SlidesView } from "./viewers/slides/SlidesView";
import { ChannelsView } from "./viewers/channels/ChannelsView";
import { SettingsView } from "./viewers/SettingsView";
import { ChannelsSettingsView } from "./viewers/settings/ChannelsSettingsView";

export function ProjectView() {
  const { setActiveNode, activeLayout } = useLayoutsStore(
    (state) => ({
      setActiveNode: state.setActiveNode,
      activeLayout: state.activeLayout,
    }),
    shallow
  );

  const handleLayoutChange = (node: MosaicNode<ViewId> | null) => {
    setActiveNode(node);
  };

  // const handleCreateNode = (): ViewId => {
  //   return "new";
  // };

  return (
    <div className={styles.container}>
      <Mosaic<ViewId>
        resize={{
          minimumPaneSizePercentage: 15,
        }}
        initialValue={activeLayout.node}
        onChange={handleLayoutChange}
        className={classNames("mosaic-blueprint-theme", Classes.DARK)}
        renderTile={(id, path) => (
          <MosaicWindow<ViewId> path={path} title={TITLE_MAP[id]}>
            {
              {
                slides: <SlidesView />,
                image: <BlendView />,
                channels: <ChannelsView />,
                settings: <ChannelsSettingsView />,
                new: <SettingsView />,
              }[id]
            }
          </MosaicWindow>
        )}
      />
    </div>
  );
}
