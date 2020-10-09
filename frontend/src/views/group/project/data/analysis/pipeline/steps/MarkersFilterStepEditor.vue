<template>
  <v-expansion-panels :value="0">
    <v-expansion-panel>
      <v-expansion-panel-header disable-icon-rotate>
        <template v-slot:actions>
          <v-tooltip bottom>
            <template v-slot:activator="{ on }">
              <v-btn icon v-on="on" @click="deleteStep(step)">
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
  @Prop(Object) step;
  @Prop(Function) deleteStep;
}
</script>
