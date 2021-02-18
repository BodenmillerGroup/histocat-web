import styles from "./ProjectView.module.scss";
import React, { useState } from "react";
import "react-mosaic-component/react-mosaic-component.css";
import { Mosaic, MosaicWindow } from "react-mosaic-component";
import classNames from "classnames";
import { Classes } from "@blueprintjs/core";

type ViewId = "a" | "b" | "c" | "new";

const TITLE_MAP: Record<ViewId, string> = {
  a: "Left Window",
  b: "Top Right Window",
  c: "Bottom Right Window",
  new: "New Window",
};

const LOCAL_STORAGE_KEY = "ProjectLayout";
const DEFAULT_LAYOUT: any = {
  direction: "row",
  first: "a",
  second: {
    direction: "column",
    first: "b",
    second: "c",
  },
};

export function ProjectView() {
  const [layout, setLayout] = useState(DEFAULT_LAYOUT);

  const saveLayout = () => {
    localStorage.setItem(LOCAL_STORAGE_KEY, JSON.stringify(layout));
  }

  const loadLayout = () => {
    const layout = localStorage.getItem(LOCAL_STORAGE_KEY);
    if (layout) {
      setLayout(JSON.parse(layout));
    } else {
      setLayout(DEFAULT_LAYOUT);
    }
  }

  const deleteLayout = () => {
    localStorage.removeItem(LOCAL_STORAGE_KEY);
  }

  return (
    <div className={styles.container}>
      <Mosaic<ViewId>
        initialValue={DEFAULT_LAYOUT}
        onChange={(node) => setLayout(node)}
        // onChange={(node) => console.log(node)}
        className={classNames("mosaic-blueprint-theme", Classes.DARK)}
        renderTile={(id, path) => (
          <MosaicWindow<ViewId> path={path} createNode={() => "new"} title={TITLE_MAP[id]}>
            <h1>{TITLE_MAP[id]}</h1>
          </MosaicWindow>
        )}
      />
    </div>
  );
}
