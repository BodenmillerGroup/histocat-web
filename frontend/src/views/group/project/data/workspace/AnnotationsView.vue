<template>
  <v-card tile>
    <v-toolbar flat dense color="grey lighten-4">
      <v-btn
        @click.stop="
          name = '';
          cellClass = '';
          addDialog = true;
        "
        :disabled="selectedCells.length === 0"
        color="primary"
        elevation="1"
        small
        >Add annotation</v-btn
      >
    </v-toolbar>
    <v-list dense class="overflow-y-auto scroll-view pa-0">
      <v-list-item-group v-model="selected" color="primary">
        <v-list-item v-for="item in annotations" :key="item[0]">
          <v-list-item-avatar size="16" :color="annotationsContext.getters.classes[item.cellClass]" />
          <v-list-item-action>
            <v-checkbox v-model="item.visible" hide-details dense />
          </v-list-item-action>
          <v-list-item-content>
            <v-list-item-title class="item">
              <span>{{ item.name }}</span>
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
                      name = item.name;
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
                  <v-btn icon small v-on="on" color="secondary lighten-2" @click.stop="deleteAnnotation(item.name)">
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
          <v-text-field label="Name" v-model="name" />
          <v-select dense :items="classes" v-model="cellClass" label="Class" />
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
          <v-text-field label="Name" v-model="name" />
          <v-select dense :items="classes" v-model="cellClass" label="Class" />
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

@Component
export default class AnnotationsView extends Vue {
  readonly annotationsContext = annotationsModule.context(this.$store);
  readonly cellsContext = cellsModule.context(this.$store);

  addDialog = false;
  editDialog = false;
  activeId: number | null = null;
  name: string | null = null;
  cellClass: string | null = null;

  selected?: any | null = null;

  get annotations() {
    return this.annotationsContext.getters.annotations;
  }

  get classes() {
    return Object.keys(this.annotationsContext.getters.classes).sort((a, b) => a.localeCompare(b));
  }

  get selectedCellIds() {
    return this.cellsContext.getters.selectedCellIds;
  }

  deleteAnnotation(name: string) {
    if (self.confirm("Do you really want to delete this annotation?")) {
      this.annotationsContext.mutations.deleteAnnotation(name);
    }
  }

  updateAnnotation() {
    this.editDialog = false;
    if (this.selected) {
      this.annotationsContext.mutations.updateAnnotation({
        prevName: this.selected.name,
        nextName: this.name!,
        cellClass: this.cellClass!,
      });
    }
  }

  addAnnotation() {
    this.addDialog = false;
    if (this.selectedCellIds) {
      this.annotationsContext.mutations.addAnnotation({
        name: this.name!,
        cellClass: this.cellClass!,
        cells: this.selectedCellIds,
      });
    }
  }
}
</script>

<style scoped>
.scroll-view {
  height: calc(33vh - 100px);
}
.item {
  display: flex;
  justify-content: space-between;
}
</style>
