<template>
  <v-data-table
    :headers="headers"
    :items="items"
    :search="search"
    hide-default-footer
    dense
    disable-pagination
    disable-sort
    no-data-text="Please create annotations first"
  >
    <template v-slot:top>
      <v-text-field v-model="search" label="Search" clearable single-line>
        <template v-slot:append>
          <v-icon>mdi-magnify</v-icon>
        </template>
      </v-text-field>
    </template>
    <template v-slot:item.threshold="props">
      <v-text-field
        v-model.number="props.item.threshold"
        type="number"
        min="0"
        max="1"
        step="0.1"
        :rules="probabilityRules"
        dense
        single-line
      />
    </template>
  </v-data-table>
</template>

<script lang="ts">
import { uniqBy } from "lodash-es";
import { Component, Vue } from "vue-property-decorator";
import { annotationsModule } from "@/modules/annotations";
import { percentFloatNumber } from "@/utils/validators";

@Component
export default class ThresholdSelector extends Vue {
  readonly annotationsContext = annotationsModule.context(this.$store);

  readonly probabilityRules = [percentFloatNumber];

  search = "";

  readonly headers = [
    {
      text: "Cell Class",
      value: "cellClass",
      align: "start",
    },
    {
      text: "Threshold",
      value: "threshold",
      filterable: false,
      align: "end",
    },
  ];

  get annotations() {
    return this.annotationsContext.getters.annotations;
  }

  get items() {
    return uniqBy(this.annotations, "cellClass").map((annotation) => {
      return {
        cellClass: annotation.cellClass,
        threshold: 0.5,
      };
    });
  }

  public getThresholds() {
    const thresholds = {};
    for (const item of this.items) {
      thresholds[item.cellClass] = item.threshold;
    }
    return thresholds;
  }
}
</script>
