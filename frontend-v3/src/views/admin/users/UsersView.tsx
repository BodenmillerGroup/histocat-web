import styles from "./UsersView.module.scss";
import { useUsersStore } from "../../../modules/users";
import shallow from "zustand/shallow";
import { useEffect, useState } from "react";
import { ActionsColumn, TextSortableColumn, NumericSortableColumn } from "../../../components/table/Core";
import { Cell, SelectionModes, Table, Utils } from "@blueprintjs/table";
import { Button, FormGroup, InputGroup } from "@blueprintjs/core";

export function UsersView() {
  const { users, getUsers } = useUsersStore((state) => ({ users: state.users, getUsers: state.getUsers }), shallow);
  const [sortedIndexMap, setSortedIndexMap] = useState<number[]>([]);
  const [filterValue, setFilterValue] = useState<string>("");

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
    window.setTimeout(() => setFilterValue((event.target as HTMLInputElement).value), 10);

  const editAction = (item: any) => {
    console.log(item);
  };

  const deleteAction = (item: any) => {
    console.log(item);
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
    new TextSortableColumn("Email", "email").getColumn(getCellData, sortColumn),
    new TextSortableColumn("Name", "name").getColumn(getCellData, sortColumn),
    new ActionsColumn("Actions").getColumn(actionsCellRenderer),
  ];

  const clearButton = <Button icon="delete" minimal={true} onClick={() => setFilterValue("")} />;

  return (
    <div className={styles.container}>
      <h2>Users</h2>
      <InputGroup
        asyncControl={true}
        leftIcon="filter"
        rightElement={filterValue ? clearButton : undefined}
        onChange={handleFilterChange}
        placeholder="Filter users..."
        value={filterValue}
        style={{ marginBottom: "20px" }}
      />
      <Table
        numRows={data.length}
        selectionModes={SelectionModes.ROWS_ONLY}
        enableRowReordering={false}
        enableColumnReordering={false}
        enableFocusedCell={false}
        enableRowResizing={false}
        defaultRowHeight={24}
        columnWidths={[50, 200, 200, 200]}
      >
        {columns}
      </Table>
    </div>
  );
}
