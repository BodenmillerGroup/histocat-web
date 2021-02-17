import styles from "./ModelsView.module.scss";
import shallow from "zustand/shallow";
import { useEffect, useState } from "react";
import { ActionsColumn, NumericSortableColumn, TextSortableColumn } from "components/table/Core";
import { Cell, SelectionModes, Table, Utils } from "@blueprintjs/table";
import { Button, InputGroup } from "@blueprintjs/core";
import { throttle } from "lodash-es";
import { EditModelDialog } from "./EditModelDialog";
import { AddModelDialog } from "./AddModelDialog";
import { useModelsStore } from "modules/models";
import { IModel } from "modules/models/models";

export function ModelsView() {
  const { ids, entities, getModels, deleteModel } = useModelsStore(
    (state) => ({
      ids: state.ids,
      entities: state.entities,
      getModels: state.getModels,
      deleteModel: state.deleteModel,
    }),
    shallow
  );
  const [sortedIndexMap, setSortedIndexMap] = useState<number[]>([]);
  const [filterValue, setFilterValue] = useState<string>("");
  const [editDialogOpen, setEditDialogOpen] = useState<boolean>(false);
  const [addDialogOpen, setAddDialogOpen] = useState<boolean>(false);
  const [activeItem, setActiveItem] = useState<IModel | null>(null);

  const models = ids.map((id) => entities[id]);
  const data = models.filter((item) => {
    const filter = filterValue.toLowerCase();
    return (
      item.name.toLowerCase().includes(filter) || (item.description && item.description.toLowerCase().includes(filter))
    );
  });

  useEffect(() => {
    getModels();
  }, [getModels]);

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

  const editAction = (item: IModel) => {
    setActiveItem(item);
    setEditDialogOpen(true);
  };

  const deleteAction = async (item: IModel) => {
    if (window.confirm("Are you sure you want to delete the model?")) {
      await deleteModel(item.id);
    }
  };

  const actionsCellRenderer = (rowIndex: number) => {
    const sortedRowIndex = sortedIndexMap[rowIndex];
    if (sortedRowIndex != null) {
      rowIndex = sortedRowIndex;
    }
    return (
      <Cell>
        <Button
          text="Edit"
          icon="edit"
          minimal={true}
          small={true}
          intent="primary"
          onClick={() => editAction(data[rowIndex])}
        />
        <Button
          text="Delete"
          icon="delete"
          minimal={true}
          small={true}
          intent="danger"
          onClick={() => deleteAction(data[rowIndex])}
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
        <h2>Models</h2>
        <Button text="Add Model" icon="add" intent="primary" onClick={() => setAddDialogOpen(true)} />
      </span>
      <InputGroup
        asyncControl={true}
        leftIcon="filter"
        rightElement={
          filterValue ? <Button icon="cross" minimal={true} onClick={() => setFilterValue("")} /> : undefined
        }
        onChange={throttle(handleFilterChange, 200, { leading: false })}
        placeholder="Filter models..."
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
        columnWidths={[50, 200, 200, 200, 200]}
      >
        {columns}
      </Table>
      <EditModelDialog model={activeItem!} isOpen={editDialogOpen} handleClose={() => setEditDialogOpen(false)} />
      <AddModelDialog isOpen={addDialogOpen} handleClose={() => setAddDialogOpen(false)} />
    </div>
  );
}
