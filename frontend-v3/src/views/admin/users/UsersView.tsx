import styles from "./UsersView.module.scss";
import { useUsersStore } from "../../../modules/users";
import shallow from "zustand/shallow";
import { useEffect, useState } from "react";
import {
  ActionsColumn,
  CheckboxSortableColumn,
  NumericSortableColumn,
  TextSortableColumn,
} from "../../../components/table/Core";
import { Cell, SelectionModes, Table, Utils } from "@blueprintjs/table";
import { Button, InputGroup } from "@blueprintjs/core";
import { EditUserDialog } from "./EditUserDialog";
import { IUserProfile } from "../../../modules/profile/models";
import { AddUserDialog } from "./AddUserDialog";

export function UsersView() {
  const { users, getUsers } = useUsersStore((state) => ({ users: state.users, getUsers: state.getUsers }), shallow);
  const [sortedIndexMap, setSortedIndexMap] = useState<number[]>([]);
  const [filterValue, setFilterValue] = useState<string>("");
  const [editDialogOpen, setEditDialogOpen] = useState<boolean>(false);
  const [addDialogOpen, setAddDialogOpen] = useState<boolean>(false);
  const [activeItem, setActiveItem] = useState<IUserProfile | null>(null);

  const data = users.filter((item) => item.name.includes(filterValue) || item.email.includes(filterValue));

  useEffect(() => {
    getUsers();
  }, [getUsers]);

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
    window.setTimeout(() => setFilterValue((event.target as HTMLInputElement).value), 100);

  const editAction = (item: IUserProfile) => {
    setActiveItem(item);
    setEditDialogOpen(true);
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
      </Cell>
    );
  };

  const columns = [
    new NumericSortableColumn("ID", "id").getColumn(getCellData, sortColumn),
    new TextSortableColumn("Email", "email").getColumn(getCellData, sortColumn),
    new TextSortableColumn("Name", "name").getColumn(getCellData, sortColumn),
    new CheckboxSortableColumn("Active", "is_active").getColumn(getCellData, sortColumn),
    new CheckboxSortableColumn("Admin", "is_admin").getColumn(getCellData, sortColumn),
    new ActionsColumn("Actions").getColumn(actionsCellRenderer),
  ];

  return (
    <div className={styles.container}>
      <span className={styles.toolbar}>
        <h2>Users</h2>
        <Button text="Add User" icon="add" intent="primary" onClick={() => setAddDialogOpen(true)} />
      </span>
      <InputGroup
        asyncControl={true}
        leftIcon="filter"
        rightElement={
          filterValue ? <Button icon="delete" minimal={true} onClick={() => setFilterValue("")} /> : undefined
        }
        onChange={handleFilterChange}
        placeholder="Filter users..."
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
        columnWidths={[50, 200, 200, 70, 70, 200]}
      >
        {columns}
      </Table>
      <EditUserDialog user={activeItem!} isOpen={editDialogOpen} handleClose={() => setEditDialogOpen(false)} />
      <AddUserDialog isOpen={addDialogOpen} handleClose={() => setAddDialogOpen(false)} />
    </div>
  );
}
