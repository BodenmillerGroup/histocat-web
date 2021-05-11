<template>
  <div>
    <table>
      <tbody>
        <tr>
          <th>CellID:</th>
          <td>{{ cell.cellId }}</td>
        </tr>
        <tr>
          <th>ObjectNumber:</th>
          <td>{{ cell.objectNumber }}</td>
        </tr>
        <tr>
          <th>AcquisitionId:</th>
          <td>{{ cell.acquisitionId }}</td>
        </tr>
      </tbody>
    </table>
  </div>
</template>

<script lang="ts">
import { Component, Prop, Vue } from "vue-property-decorator";
import { ICell } from "@/modules/cells/models";

@Component
export default class Tooltip extends Vue {
  @Prop(Object) readonly cell!: ICell;
}
</script>
