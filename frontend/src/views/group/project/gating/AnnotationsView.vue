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
      <v-btn
        @click.stop="predictDialog = true"
        :disabled="!activeDataset || annotations.length === 0"
        color="primary"
        elevation="1"
        x-small
        class="ml-2"
        >Predict</v-btn
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
    <v-dialog v-model="predictDialog" scrollable max-width="900px">
      <v-card>
        <v-card-title>Predict cell classes</v-card-title>
        <v-card-text>
          <v-row>
            <v-col cols="8"><ChannelSelector ref="channelSelector" /></v-col>
            <v-col cols="4">
              <ThresholdSelector ref="thresholdSelector" />
              <v-text-field
                label="n Estimators"
                v-model.number="nEstimators"
                type="number"
                min="1"
                step="1"
                dense
                :rules="nEstimatorsRules"
                class="ma-8"
              />
            </v-col>
          </v-row>
        </v-card-text>
        <v-card-actions>
          <v-spacer />
          <v-btn text @click="predictDialog = false">Cancel</v-btn>
          <v-btn color="primary" text @click="trainCellClassifier">Predict</v-btn>
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
import { analysisModule } from "@/modules/analysis";
import { datasetsModule } from "@/modules/datasets";
import ChannelSelector from "@/views/group/project/gating/ChannelSelector.vue";
import ThresholdSelector from "@/views/group/project/gating/ThresholdSelector.vue";
import { required, positiveNumber } from "@/utils/validators";

@Component({
  components: { ThresholdSelector, ChannelSelector },
})
export default class AnnotationsView extends Vue {
  readonly annotationsContext = annotationsModule.context(this.$store);
  readonly analysisContext = analysisModule.context(this.$store);
  readonly cellsContext = cellsModule.context(this.$store);
  readonly datasetsContext = datasetsModule.context(this.$store);

  readonly nEstimatorsRules = [required, positiveNumber];

  addDialog = false;
  editDialog = false;
  predictDialog = false;
  activeId: number | null = null;
  cellClass: string | null = null;

  nEstimators = 100;

  selected?: any | null = null;

  get activeDataset() {
    return this.datasetsContext.getters.activeDataset;
  }

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
        cellIds: this.selectedCellIds,
      });
    }
  }

  trainCellClassifier() {
    this.predictDialog = false;
    const channels = (this.$refs.channelSelector as any).getChannels();
    const thresholds = (this.$refs.thresholdSelector as any).getThresholds();
    this.analysisContext.actions.classifyCells({
      channels: channels,
      thresholds: thresholds,
      nEstimators: this.nEstimators,
    });
  }
}
</script>

<style scoped>
.scroll-view {
  height: calc(33vh - 100px);
}
</style>
