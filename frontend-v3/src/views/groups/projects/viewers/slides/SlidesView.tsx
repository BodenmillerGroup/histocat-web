import styles from "./SlidesView.module.scss";
import shallow from "zustand/shallow";
import React, { useEffect, useState } from "react";
import { Alert, Button, Classes, InputGroup, Intent, ITreeNode, Menu, MenuItem, Tree } from "@blueprintjs/core";
import { throttle } from "lodash-es";
import { AddSlideDialog } from "./AddSlideDialog";
import { ContextMenu2, Popover2 } from "@blueprintjs/popover2";
import { useProjectsStore } from "modules/projects";
import { IAcquisition, ISlide } from "modules/projects/models";
import { useDatasetsStore } from "modules/datasets";
import ReactJson from "react-json-view";

function forEachNode(nodes: ITreeNode[], callback: (node: ITreeNode) => void) {
  if (nodes == null) {
    return;
  }

  for (const node of nodes) {
    callback(node);
    forEachNode(node.childNodes!, callback);
  }
}

export function SlidesView() {
  const { projectData, setActiveAcquisitionId, deleteSlide } = useProjectsStore(
    (state) => ({
      projectData: state.projectData,
      setActiveAcquisitionId: state.setActiveAcquisitionId,
      deleteSlide: state.deleteSlide,
    }),
    shallow
  );
  const { activeDataset } = useDatasetsStore(
    (state) => ({
      activeDataset: state.getActiveDataset(),
    }),
    shallow
  );
  const [filterValue, setFilterValue] = useState<string>("");
  const [alertOpen, setAlertOpen] = useState<boolean>(false);
  const [addDialogOpen, setAddDialogOpen] = useState<boolean>(false);
  const [activeItem, setActiveItem] = useState<ISlide | null>(null);
  let data: ITreeNode[] = projectData!.slides!.map((slide) => {
        const acquisitions = slide.acquisitions.map((acquisition) => {
          let hasMask = false;
          if (activeDataset && activeDataset.meta.masks) {
            hasMask = !!activeDataset.meta.masks[acquisition.id];
          }
          return {
            hasCaret: false,
            label: acquisition.description,
            id: `acquisition-${acquisition.id}`,
            disabled: !acquisition.is_valid,
            icon: hasMask && "heatmap",
            nodeData: acquisition,
            secondaryLabel: (
              <Popover2
                content={
                  <ReactJson src={acquisition.meta} displayDataTypes={false} displayObjectSize={false} theme="apathy" />
                }
                placement="right"
              >
                <Button minimal={true} small={true} icon="info-sign" />
              </Popover2>
            ),
          };
        });
        const panoramas = slide.panoramas.map((panorama) => {
          return {
            hasCaret: false,
            label: panorama.description,
            id: `panorama-${panorama.id}`,
            disabled: panorama.image_type === "Default",
            nodeData: panorama,
            secondaryLabel: (
              <Popover2
                content={
                  <ReactJson src={panorama.meta} displayDataTypes={false} displayObjectSize={false} theme="apathy" />
                }
                placement="right"
              >
                <Button minimal={true} small={true} icon="info-sign" />
              </Popover2>
            ),
          };
        });
        const panoramasRoot = {
          label: "Panorama Images",
          id: `slide-${slide.id}-panoramas`,
          childNodes: panoramas,
          isExpanded: false,
        };
        return {
          icon: "folder-close",
          label: (
            <ContextMenu2
              content={
                <Menu>
                  <MenuItem text="Delete slide..." icon="delete" intent="danger" onClick={() => setAlertOpen(true)} />
                </Menu>
              }
            >
              {slide.name}
            </ContextMenu2>
          ),
          childNodes: [panoramasRoot as any].concat(acquisitions),
          isExpanded: true,
          id: `slide-${slide.id}`,
          nodeData: slide,
          secondaryLabel: (
            <Popover2
              content={<ReactJson src={slide.meta} displayDataTypes={false} displayObjectSize={false} theme="apathy" />}
              placement="right"
            >
              <Button minimal={true} small={true} icon="info-sign" />
            </Popover2>
          ),
        };
      });
  data = data.filter((item) => {
    const filter = filterValue.toLowerCase();
    return (item.nodeData as any).name.toLowerCase().includes(filter);
  });
  const [nodes, setNodes] = useState<ITreeNode[]>(data);

  useEffect(() => {
    if (data.length !== nodes.length) {
      setNodes([...(data as ITreeNode[])]);
    }
  }, [data, nodes]);

  const handleNodeClick = (node: ITreeNode, _nodePath: number[], e: React.MouseEvent<HTMLElement>) => {
    const originallySelected = node.isSelected;
    if (!e.shiftKey) {
      forEachNode(nodes, (n) => (n.isSelected = false));
    }
    node.isSelected = originallySelected == null ? true : !originallySelected;
    setNodes([...nodes]);
    if (node.nodeData && node.nodeData.hasOwnProperty("channels")) {
      setActiveAcquisitionId(node.isSelected ? (node.nodeData as IAcquisition).id : null);
    }
  };

  const handleNodeCollapse = (nodeData: ITreeNode) => {
    nodeData.isExpanded = false;
    setNodes([...nodes]);
  };

  const handleNodeExpand = (nodeData: ITreeNode) => {
    nodeData.isExpanded = true;
    setNodes([...nodes]);
  };

  const handleFilterChange = (event: React.FormEvent<HTMLElement>) =>
    setFilterValue((event.target as HTMLInputElement).value);

  const deleteAction = async () => {
    if (activeItem) {
      await deleteSlide(activeItem.id);
    }
    setAlertOpen(false);
  };

  return (
    <div className={styles.container}>
      <span className={styles.toolbar}>
        <InputGroup
          asyncControl={true}
          leftIcon="filter"
          rightElement={
            filterValue ? <Button icon="cross" minimal={true} onClick={() => setFilterValue("")} /> : undefined
          }
          onChange={throttle(handleFilterChange, 200, { leading: false })}
          placeholder="Filter slides..."
          value={filterValue}
          className={styles.filter}
        />
        <Button text="Upload" icon="upload" intent={Intent.PRIMARY} onClick={() => setAddDialogOpen(true)} />
      </span>
      <Tree
        contents={nodes}
        onNodeClick={handleNodeClick}
        onNodeCollapse={handleNodeCollapse}
        onNodeExpand={handleNodeExpand}
        onNodeContextMenu={(node) => setActiveItem(node.nodeData as ISlide)}
        className={styles.scrollable}
      />
      <AddSlideDialog isOpen={addDialogOpen} handleClose={() => setAddDialogOpen(false)} />
      <Alert
        className={Classes.DARK}
        cancelButtonText="Cancel"
        canEscapeKeyCancel={true}
        confirmButtonText="Delete"
        icon="trash"
        intent={Intent.DANGER}
        isOpen={alertOpen}
        onCancel={() => setAlertOpen(false)}
        onConfirm={deleteAction}
      >
        <p>Are you sure you want to delete the slide?</p>
      </Alert>
    </div>
  );
}
