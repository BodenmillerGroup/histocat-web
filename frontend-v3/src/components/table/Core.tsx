import { Cell, Column, ColumnHeaderCell, ICellRenderer } from "@blueprintjs/table";
import { Menu, MenuItem } from "@blueprintjs/core";

export type ICellLookup = (rowIndex: number, accesor: string) => any;
export type ISortCallback = (accessor: string, comparator: (a: any, b: any) => number) => void;

export class TextSortableColumn {
  constructor(protected name: string, protected accessor: string) {}

  private renderMenu(sortColumn: ISortCallback) {
    const sortAsc = () => sortColumn(this.accessor, (a, b) => this.compare(a, b));
    const sortDesc = () => sortColumn(this.accessor, (a, b) => this.compare(b, a));
    return (
      <Menu>
        <MenuItem icon="sort-asc" onClick={sortAsc} text="Sort Asc" />
        <MenuItem icon="sort-desc" onClick={sortDesc} text="Sort Desc" />
      </Menu>
    );
  }

  private compare(a: string, b: string) {
    return a.toString().localeCompare(b);
  }

  public getColumn(getCellData: ICellLookup, sortColumn: ISortCallback) {
    const cellRenderer = (rowIndex: number, columnIndex: number) => <Cell>{getCellData(rowIndex, this.accessor)}</Cell>;
    const menuRenderer = this.renderMenu.bind(this, sortColumn);
    const columnHeaderCellRenderer = () => <ColumnHeaderCell name={this.name} menuRenderer={menuRenderer} />;
    return (
      <Column
        cellRenderer={cellRenderer}
        columnHeaderCellRenderer={columnHeaderCellRenderer}
        key={this.name}
        name={this.name}
      />
    );
  }
}

export class NumericSortableColumn {
  constructor(protected name: string, protected accessor: string) {}

  private renderMenu(sortColumn: ISortCallback) {
    const sortAsc = () => sortColumn(this.accessor, (a, b) => this.compare(a, b));
    const sortDesc = () => sortColumn(this.accessor, (a, b) => this.compare(b, a));
    return (
      <Menu>
        <MenuItem icon="sort-asc" onClick={sortAsc} text="Sort Asc" />
        <MenuItem icon="sort-desc" onClick={sortDesc} text="Sort Desc" />
      </Menu>
    );
  }

  private compare(a: number, b: number) {
    return a - b;
  }

  public getColumn(getCellData: ICellLookup, sortColumn: ISortCallback) {
    const cellRenderer = (rowIndex: number, columnIndex: number) => <Cell>{getCellData(rowIndex, this.accessor)}</Cell>;
    const menuRenderer = this.renderMenu.bind(this, sortColumn);
    const columnHeaderCellRenderer = () => <ColumnHeaderCell name={this.name} menuRenderer={menuRenderer} />;
    return (
      <Column
        cellRenderer={cellRenderer}
        columnHeaderCellRenderer={columnHeaderCellRenderer}
        key={this.name}
        name={this.name}
      />
    );
  }
}

export class ActionsColumn {
  constructor(protected name: string) {}

  public getColumn(cellRenderer: ICellRenderer) {
    const columnHeaderCellRenderer = () => <ColumnHeaderCell name={this.name} />;
    return (
      <Column
        cellRenderer={cellRenderer}
        columnHeaderCellRenderer={columnHeaderCellRenderer}
        key={this.name}
        name={this.name}
      />
    );
  }
}
