<template>
  <v-card tile>
    <v-toolbar flat dense color="grey lighten-4">
      <v-btn
        @click.stop="
          name = '';
          cellClass = '';
          addDialog = true;
        "
        :disabled="selectedCellIds.length === 0"
        color="primary"
        elevation="1"
        x-small
        >Add annotation</v-btn
      >
    </v-toolbar>
    <v-list dense class="overflow-y-auto scroll-view pa-0">
      <v-list-item-group v-model="selected" color="primary">
        <v-list-item v-for="(item, index) in annotations" :key="index">
          <v-list-item-avatar size="16" :color="cellClasses[item.cellClass]" />
          <v-list-item-action>
            <v-checkbox v-model="item.visible" hide-details dense />
          </v-list-item-action>
          <v-list-item-content>
            <v-list-item-title>
              <span>{{ item.cellClass }}</span>
            </v-list-item-title>
          </v-list-item-content>
          <v-list-item-action>
            <v-row>
              <v-tooltip bottom>
                <template v-slot:activator="{ on }">
                  <v-btn
                    icon
                    small
                    v-on="on"
                    color="primary lighten-2"
                    @click.stop="
                      selected = item;
                      cellClass = item.cellClass;
                      editDialog = true;
                    "
                  >
                    <v-icon small>mdi-pencil-outline</v-icon>
                  </v-btn>
                </template>
                <span>Edit annotation</span>
              </v-tooltip>
              <v-tooltip bottom>
                <template v-slot:activator="{ on }">
                  <v-btn icon small v-on="on" color="secondary lighten-2" @click.stop="deleteAnnotation(item)">
                    <v-icon small>mdi-delete-outline</v-icon>
                  </v-btn>
                </template>
                <span>Delete annotation</span>
              </v-tooltip>
            </v-row>
          </v-list-item-action>
        </v-list-item>
      </v-list-item-group>
    </v-list>
    <v-dialog v-model="addDialog" scrollable max-width="300px">
      <v-card>
        <v-card-title>Add annotation</v-card-title>
        <v-card-text>
          <v-select :items="cellClassItems" v-model="cellClass" label="Class" />
        </v-card-text>
        <v-card-actions>
          <v-spacer />
          <v-btn text @click="addDialog = false">Cancel</v-btn>
          <v-btn color="primary" text @click="addAnnotation">Save</v-btn>
        </v-card-actions>
      </v-card>
    </v-dialog>
    <v-dialog v-model="editDialog" scrollable max-width="300px">
      <v-card>
        <v-card-title>Edit annotation</v-card-title>
        <v-card-text>
          <v-select :items="cellClassItems" v-model="cellClass" label="Class" />
        </v-card-text>
        <v-card-actions>
          <v-spacer />
          <v-btn text @click="editDialog = false">Cancel</v-btn>
          <v-btn color="primary" text @click="updateAnnotation">Save</v-btn>
        </v-card-actions>
      </v-card>
    </v-dialog>
  </v-card>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import { annotationsModule } from "@/modules/annotations";
import { cellsModule } from "@/modules/cells";
import { IAnnotation } from "@/modules/annotations/models";

@Component
export default class AnnotationsView extends Vue {
  readonly annotationsContext = annotationsModule.context(this.$store);
  readonly cellsContext = cellsModule.context(this.$store);

  addDialog = false;
  editDialog = false;
  activeId: number | null = null;
  cellClass: string | null = null;

  selected?: any | null = null;

  get annotations() {
    return this.annotationsContext.getters.annotations;
  }

  get cellClasses() {
    return this.annotationsContext.getters.cellClasses;
  }

  get cellClassItems() {
    return Object.keys(this.cellClasses).sort((a, b) => a.localeCompare(b));
  }

  get selectedCellIds() {
    return this.cellsContext.getters.selectedCellIds;
  }

  deleteAnnotation(item: IAnnotation) {
    if (self.confirm("Do you really want to delete this annotation?")) {
      const index = this.annotations.indexOf(item);
      this.annotationsContext.mutations.deleteAnnotation(index);
    }
  }

  updateAnnotation() {
    this.editDialog = false;
    if (this.selected) {
      const index = this.annotations.indexOf(this.selected);
      this.annotationsContext.mutations.updateAnnotation({
        index: index,
        cellClass: this.cellClass!,
      });
    }
  }

  addAnnotation() {
    this.addDialog = false;
    if (this.selectedCellIds) {
      this.annotationsContext.mutations.addAnnotation({
        cellClass: this.cellClass!,
        cells: this.selectedCellIds,
      });
    }
    const colors = this.selectedCellIds.map(cellId => this.cellClasses[this.cellClass!]);
    console.log(colors)
    this.cellsContext.mutations.updateCellsByColors({
      cellIds: this.selectedCellIds,
      colors: {
        type: "clustering",
        name: "Annotations",
        data: colors,
      }
    });
  }
}
</script>

<style scoped>
.scroll-view {
  height: calc(33vh - 100px);
}
</style>
