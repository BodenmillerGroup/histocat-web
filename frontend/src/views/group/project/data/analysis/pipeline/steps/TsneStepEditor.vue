<template>
  <v-expansion-panels :value="0">
    <v-expansion-panel>
      <v-expansion-panel-header>t-SNE</v-expansion-panel-header>
      <v-expansion-panel-content>
        <v-card flat>
          <v-card-text>
            Lorem ipsum dolor sit amet, no nam oblique veritus. Commune scaevola imperdiet nec ut, sed euismod convenire
            principes at. Est et nobis iisque percipit, an vim zril disputando voluptatibus, vix an salutandi
            sententiae.
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
import {Component, Prop, Vue} from "vue-property-decorator";

@Component
export default class TsneStepEditor extends Vue {
  @Prop(Object) step;
  @Prop(Function) deleteStep;
}
</script>
