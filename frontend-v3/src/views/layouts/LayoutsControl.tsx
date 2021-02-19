import { Button, ControlGroup, HTMLSelect } from "@blueprintjs/core";
import shallow from "zustand/shallow";
import { useLayoutsStore } from "modules/layouts";
import { ChangeEvent, useState } from "react";
import React from "react";
import { AddLayoutDialog } from "./AddLayoutDialog";

export function LayoutsControl() {
  const { layouts, loadLayout, deleteLayout, activeLayout } = useLayoutsStore(
    (state) => ({
      layouts: state.layouts,
      loadLayout: state.loadLayout,
      deleteLayout: state.deleteLayout,
      activeLayout: state.activeLayout,
    }),
    shallow
  );
  const [addDialogOpen, setAddDialogOpen] = useState<boolean>(false);

  const handleDeleteLayout = () => {
    deleteLayout(activeLayout.name);
  };

  const handleOnChange = (event: ChangeEvent<HTMLSelectElement>) => {
    loadLayout(event.currentTarget.value);
  };

  return (
    <React.Fragment>
      <ControlGroup vertical={false}>
        <HTMLSelect value={activeLayout.name} onChange={handleOnChange}>
          {layouts.map((layout) => {
            return <option>{layout.name}</option>;
          })}
        </HTMLSelect>
        <Button icon="add" onClick={() => setAddDialogOpen(true)} />
        <Button icon="delete" disabled={activeLayout.isDefault} onClick={handleDeleteLayout} />
      </ControlGroup>
      <AddLayoutDialog isOpen={addDialogOpen} handleClose={() => setAddDialogOpen(false)} />
    </React.Fragment>
  );
}
