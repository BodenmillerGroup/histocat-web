<template>
  <div class="cell-classes-view">
    <v-toolbar flat dense color="grey lighten-4">
      <v-btn
        @click.stop="
          name = '';
          color = '';
          addDialog = true;
        "
        color="primary"
        elevation="1"
        x-small
        tile
        >Add cell class</v-btn
      >
      <v-btn @click="resetCellClasses" color="primary" elevation="1" x-small tile class="ml-2">Reset</v-btn>
    </v-toolbar>
    <v-list dense class="pa-0">
      <v-list-item v-for="item in cellClasses" :key="item[0]">
        <v-list-item-avatar size="16" :color="item[1]" />
        <v-list-item-content>
          <v-list-item-title>{{ item[0] }}</v-list-item-title>
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
                    prevName = item[0];
                    name = item[0];
                    color = item[1];
                    editDialog = true;
                  "
                >
                  <v-icon small>mdi-pencil-outline</v-icon>
                </v-btn>
              </template>
              <span>Edit cell class</span>
            </v-tooltip>
            <v-tooltip bottom>
              <template v-slot:activator="{ on }">
                <v-btn icon small v-on="on" color="secondary lighten-2" @click.stop="deleteCellClass(item[0])">
                  <v-icon small>mdi-delete-outline</v-icon>
                </v-btn>
              </template>
              <span>Delete cell class</span>
            </v-tooltip>
          </v-row>
        </v-list-item-action>
      </v-list-item>
    </v-list>
    <v-dialog v-model="addDialog" scrollable max-width="300px">
      <v-card>
        <v-card-title>Add cell class</v-card-title>
        <v-card-text>
          <v-text-field label="Name" v-model="name" />
          <v-color-picker v-model="color"></v-color-picker>
        </v-card-text>
        <v-card-actions>
          <v-spacer />
          <v-btn text @click="addDialog = false">Cancel</v-btn>
          <v-btn color="primary" text @click="addCellClass">Save</v-btn>
        </v-card-actions>
      </v-card>
    </v-dialog>
    <v-dialog v-model="editDialog" scrollable max-width="300px">
      <v-card>
        <v-card-title>Edit cell class</v-card-title>
        <v-card-text>
          <v-text-field label="Name" v-model="name" />
          <v-color-picker v-model="color"></v-color-picker>
        </v-card-text>
        <v-card-actions>
          <v-spacer />
          <v-btn text @click="editDialog = false">Cancel</v-btn>
          <v-btn color="primary" text @click="updateCellClass">Save</v-btn>
        </v-card-actions>
      </v-card>
    </v-dialog>
  </div>
</template>

<script lang="ts">
import { Component, Vue } from "vue-property-decorator";
import { annotationsModule } from "@/modules/annotations";
import { uiModule } from "@/modules/ui";
import { cellsModule } from "@/modules/cells";
import { projectsModule } from "@/modules/projects";

@Component
export default class CellClassesView extends Vue {
  readonly uiContext = uiModule.context(this.$store);
  readonly annotationsContext = annotationsModule.context(this.$store);
  readonly cellsContext = cellsModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);

  addDialog = false;
  editDialog = false;
  activeId: number | null = null;
  name: string | null = null;
  prevName: string | null = null;
  color: string | null = null;

  get cellClasses() {
    return Object.entries(this.annotationsContext.getters.cellClasses).sort((a, b) => a[0].localeCompare(b[0]));
  }

  deleteCellClass(name: string) {
    if (self.confirm("Do you really want to delete this cell class?")) {
      this.annotationsContext.mutations.deleteCellClass(name);
    }
  }

  addCellClass() {
    this.addDialog = false;
    this.annotationsContext.mutations.addCellClass({
      name: this.name!,
      color: this.color!,
    });
  }

  updateCellClass() {
    this.editDialog = false;
    if (this.prevName) {
      this.annotationsContext.mutations.updateCellClass({
        prevName: this.prevName,
        nextName: this.name!,
        color: this.color!,
      });

      if (this.cellsContext.getters.heatmap && this.cellsContext.getters.heatmap.type === "annotation") {
        this.projectsContext.actions.getAnnotationData();
        if (this.uiContext.getters.showMask) {
          this.projectsContext.actions.getChannelStackImage();
        }
      }
    }
  }

  resetCellClasses() {
    this.annotationsContext.mutations.resetCellClasses();
  }
}
</script>

<style scoped>
.cell-classes-view {
  width: 100%;
  height: 100%;
  overflow-y: auto;
}
</style>
