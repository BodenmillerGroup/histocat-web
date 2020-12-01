<template>
  <v-expansion-panels :value="0">
    <v-expansion-panel>
      <v-expansion-panel-header disable-icon-rotate>
        <template v-slot:actions>
          <v-tooltip bottom>
            <template v-slot:activator="{ on }">
              <v-btn icon v-on="on" @click.stop="moveUpStep(step)">
                <v-icon color="primary">mdi-arrow-up-bold-outline</v-icon>
              </v-btn>
            </template>
            <span>Move step up</span>
          </v-tooltip>
          <v-tooltip bottom>
            <template v-slot:activator="{ on }">
              <v-btn icon v-on="on" @click.stop="moveDownStep(step)">
                <v-icon color="primary">mdi-arrow-down-bold-outline</v-icon>
              </v-btn>
            </template>
            <span>Move step down</span>
          </v-tooltip>
          <v-tooltip bottom>
            <template v-slot:activator="{ on }">
              <v-btn icon v-on="on" @click.stop="deleteStep(step)">
                <v-icon color="secondary">mdi-delete-outline</v-icon>
              </v-btn>
            </template>
            <span>Delete step</span>
          </v-tooltip>
        </template>
        Markers Filter
      </v-expansion-panel-header>
      <v-expansion-panel-content>
        <v-card flat>
          <MarkersSelector :step="step" />
        </v-card>
      </v-expansion-panel-content>
    </v-expansion-panel>
  </v-expansion-panels>
</template>

<script lang="ts">
import { Component, Prop, Vue } from "vue-property-decorator";
import MarkersSelector from "@/views/group/project/data/analysis/pipeline/steps/MarkersSelector.vue";

@Component({
  components: { MarkersSelector },
})
export default class MarkersFilterStepEditor extends Vue {
  @Prop({ type: Object, required: true }) step;
  @Prop({ type: Function, required: true }) deleteStep;
  @Prop({ type: Function, required: true }) moveUpStep;
  @Prop({ type: Function, required: true }) moveDownStep;
}
</script>
