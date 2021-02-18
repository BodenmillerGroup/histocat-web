import styles from "./ProjectsView.module.scss";
import shallow from "zustand/shallow";
import React, { useState } from "react";
import { ActionsColumn, NumericSortableColumn, TextSortableColumn } from "components/table/Core";
import { Cell, SelectionModes, Table, Utils } from "@blueprintjs/table";
import { Alert, Button, Classes, InputGroup, Intent } from "@blueprintjs/core";
import { throttle } from "lodash-es";
import { useProjectsStore } from "modules/projects";
import { IProject } from "modules/projects/models";
import { EditProjectDialog } from "./EditProjectDialog";
import { AddProjectDialog } from "./AddProjectDialog";
import { useHistory, useRouteMatch } from "react-router-dom";

export function ProjectsView() {
  const history = useHistory();
  const { url } = useRouteMatch();
  const { ids, entities, projectsTags, deleteProject } = useProjectsStore(
    (state) => ({
      ids: state.ids,
      entities: state.entities,
      projectsTags: state.projectsTags,
      deleteProject: state.deleteProject,
    }),
    shallow
  );
  const [alertOpen, setAlertOpen] = useState<boolean>(false);
  const [sortedIndexMap, setSortedIndexMap] = useState<number[]>([]);
  const [filterValue, setFilterValue] = useState<string>("");
  const [editDialogOpen, setEditDialogOpen] = useState<boolean>(false);
  const [addDialogOpen, setAddDialogOpen] = useState<boolean>(false);
  const [activeItem, setActiveItem] = useState<IProject | null>(null);

  const projects = ids.map((id) => entities[id]);
  const data = projects.filter((item) => {
    const filter = filterValue.toLowerCase();
    return (
      item.name.toLowerCase().includes(filter) || (item.description && item.description.toLowerCase().includes(filter))
    );
  });

  const getCellData = (rowIndex: number, accessor: string) => {
    const sortedRowIndex = sortedIndexMap[rowIndex];
    if (sortedRowIndex != null) {
      rowIndex = sortedRowIndex;
    }
    return (data[rowIndex] as any)[accessor];
  };

  const sortColumn = (accessor: string, comparator: (a: any, b: any) => number) => {
    const sortedIndexMap = Utils.times(data.length, (i: number) => i);
    sortedIndexMap.sort((a: number, b: number) => {
      return comparator((data[a] as any)[accessor], (data[b] as any)[accessor]);
    });
    setSortedIndexMap(sortedIndexMap);
  };

  const handleFilterChange = (event: React.FormEvent<HTMLElement>) =>
    setFilterValue((event.target as HTMLInputElement).value);

  const editAction = (item: IProject) => {
    setActiveItem(item);
    setEditDialogOpen(true);
  };

  const deleteAction = async () => {
    if (activeItem) {
      await deleteProject(activeItem.id);
    }
    setAlertOpen(false);
  };

  const actionsCellRenderer = (rowIndex: number) => {
    const sortedRowIndex = sortedIndexMap[rowIndex];
    if (sortedRowIndex != null) {
      rowIndex = sortedRowIndex;
    }
    return (
      <Cell>
        <Button
          text="Open"
          icon="document-open"
          minimal={true}
          small={true}
          intent={Intent.PRIMARY}
          onClick={() => history.push(`${url}/projects/${data[rowIndex].id}`)}
        />
        <Button text="Edit" icon="edit" minimal={true} small={true} onClick={() => editAction(data[rowIndex])} />
        <Button
          text="Delete"
          icon="delete"
          minimal={true}
          small={true}
          intent={Intent.DANGER}
          onClick={() => {
            setActiveItem(data[rowIndex]);
            setAlertOpen(true);
          }}
        />
      </Cell>
    );
  };

  const columns = [
    new NumericSortableColumn("ID", "id").getColumn(getCellData, sortColumn),
    new TextSortableColumn("Name", "name").getColumn(getCellData, sortColumn),
    new TextSortableColumn("Description", "description").getColumn(getCellData, sortColumn),
    new TextSortableColumn("Date", "created_at").getColumn(getCellData, sortColumn),
    new ActionsColumn("Actions").getColumn(actionsCellRenderer),
  ];

  return (
    <div className={styles.container}>
      <span className={styles.toolbar}>
        <h2>Projects</h2>
        <Button text="Add Project" icon="add" intent="primary" onClick={() => setAddDialogOpen(true)} />
      </span>
      <InputGroup
        asyncControl={true}
        leftIcon="filter"
        rightElement={
          filterValue ? <Button icon="cross" minimal={true} onClick={() => setFilterValue("")} /> : undefined
        }
        onChange={throttle(handleFilterChange, 200, { leading: false })}
        placeholder="Filter projects..."
        value={filterValue}
        className={styles.filter}
      />
      <Table
        numRows={data.length}
        selectionModes={SelectionModes.ROWS_ONLY}
        enableRowReordering={false}
        enableColumnReordering={false}
        enableFocusedCell={false}
        enableRowResizing={false}
        defaultRowHeight={24}
        columnWidths={[50, 200, 300, 200, 235]}
      >
        {columns}
      </Table>
      {activeItem && (
        <EditProjectDialog
          project={activeItem}
          projectsTags={projectsTags}
          isOpen={editDialogOpen}
          handleClose={() => {
            setEditDialogOpen(false);
            setActiveItem(null);
          }}
        />
      )}
      <AddProjectDialog
        projectsTags={projectsTags}
        isOpen={addDialogOpen}
        handleClose={() => setAddDialogOpen(false)}
      />
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
        <p>Are you sure you want to delete the project?</p>
      </Alert>
    </div>
  );
}
