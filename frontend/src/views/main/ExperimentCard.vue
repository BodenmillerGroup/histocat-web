<template>
  <v-card tile class="ma-6 pa-1">
    <v-card-title>
      <v-row no-gutters>
        <v-col>
          <h5 class="headline">{{ experiment.name }}</h5>
          <span class="caption"><v-icon small>mdi-calendar-outline</v-icon> {{ createdAt }}</span>
        </v-col>
      </v-row>
    </v-card-title>
    <v-card-text v-if="experiment.description">
      {{ experiment.description }}
    </v-card-text>
    <v-card-text v-if="experiment.tags">
      <v-chip :key="item" v-for="item in experiment.tags" label small class="mr-1">
        <v-icon small left>mdi-tag-outline</v-icon>
        {{ item }}
      </v-chip>
    </v-card-text>
    <v-card-actions>
      <v-btn color="primary" :to="{ name: 'main-experiment', params: { experimentId: experiment.id } }">
        Open
      </v-btn>
      <v-spacer></v-spacer>
      <v-tooltip bottom v-if="isOwner || hasAdminAccess">
        <template v-slot:activator="{ on }">
          <v-btn icon v-on="on" :to="{ name: 'main-experiment-edit', params: { experimentId: experiment.id } }">
            <v-icon>mdi-pencil</v-icon>
          </v-btn>
        </template>
        <span>Edit experiment</span>
      </v-tooltip>
      <v-tooltip bottom v-if="isOwner">
        <template v-slot:activator="{ on }">
          <v-btn icon v-on="on" :to="{ name: 'main-experiment-share', params: { experimentId: experiment.id } }">
            <v-icon>mdi-share-variant</v-icon>
          </v-btn>
        </template>
        <span>Share experiment</span>
      </v-tooltip>
      <v-tooltip bottom v-if="isOwner || hasAdminAccess">
        <template v-slot:activator="{ on }">
          <v-btn icon color="secondary" v-on="on" @click="deleteExperiment($event, experiment.id)">
            <v-icon>mdi-delete-outline</v-icon>
          </v-btn>
        </template>
        <span>Delete experiment</span>
      </v-tooltip>
    </v-card-actions>
  </v-card>
</template>

<script lang="ts">
import { experimentModule } from "@/modules/experiment";
import { IExperiment } from "@/modules/experiment/models";
import { IUserProfile } from "@/modules/user/models";
import { Component, Prop, Vue } from "vue-property-decorator";
import { mainModule } from "@/modules/main";

@Component
export default class ExperimentCard extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly experimentContext = experimentModule.context(this.$store);

  @Prop(Object) experiment!: IExperiment;
  @Prop(Object) user!: IUserProfile;

  get createdAt() {
    return new Date(this.experiment.created_at).toUTCString();
  }

  get isOwner() {
    return this.experiment.user_id === this.user.id;
  }

  get hasAdminAccess() {
    return this.mainContext.getters.hasAdminAccess;
  }

  async deleteExperiment(event, id: number) {
    if (self.confirm("Are you sure you want to delete the experiment?")) {
      await this.experimentContext.actions.deleteExperiment(id);
    }
  }
}
</script>
