<template>
  <v-expansion-panels :value="0">
    <v-expansion-panel>
      <v-expansion-panel-header>Transformation</v-expansion-panel-header>
      <v-expansion-panel-content>
        <v-card flat>
          <v-card-text>
            <v-text-field
              label="Co-factor"
              v-model.number="step.cofactor"
              type="number"
              :min="0"
              :step="1"
              :rules="cofactorRules"
              class="number-input"
            />
          </v-card-text>
          <v-card-actions>
            <v-tooltip bottom>
              <template v-slot:activator="{ on }">
                <v-btn icon color="secondary" v-on="on" @click="deleteStep(step)">
                  <v-icon>mdi-delete-outline</v-icon>
                </v-btn>
              </template>
              <span>Delete step</span>
            </v-tooltip>
          </v-card-actions>
        </v-card>
      </v-expansion-panel-content>
    </v-expansion-panel>
  </v-expansion-panels>
</template>

<script lang="ts">
import { Component, Prop, Vue } from "vue-property-decorator";
import { required } from "@/utils/validators";

@Component
export default class TransformationStepEditor extends Vue {
  @Prop(Object) step;
  @Prop(Function) deleteStep;

  readonly cofactorRules = [required];
}
</script>

<style scoped>
.number-input {
  width: 50px;
}
</style>
