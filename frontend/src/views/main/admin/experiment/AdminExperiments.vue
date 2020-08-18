<template>
  <div>
    <v-toolbar dense light>
      <v-toolbar-title>
        Manage Experiments
      </v-toolbar-title>
      <v-spacer></v-spacer>
      <v-btn small color="primary" to="/main/experiments/create">Create Experiment</v-btn>
    </v-toolbar>

    <v-data-table :headers="headers" :items="experiments">
      <template v-slot:item.action="{ item }">
        <v-tooltip bottom>
          <template v-slot:activator="{ on }">
            <v-btn v-on="on" icon :to="{ name: 'main-admin-experiments-edit', params: { experimentId: item.id } }">
              <v-icon>mdi-pencil</v-icon>
            </v-btn>
          </template>
          <span>Edit</span>
        </v-tooltip>
        <v-tooltip bottom>
          <template v-slot:activator="{ on }">
            <v-btn v-on="on" icon @click="deleteExperiment($event, item.id)">
              <v-icon>mdi-delete</v-icon>
            </v-btn>
          </template>
          <span>Delete</span>
        </v-tooltip>
      </template>
    </v-data-table>
  </div>
</template>

<script lang="ts">
import { experimentModule } from "@/modules/experiment";
import { Component, Vue } from "vue-property-decorator";

@Component
export default class AdminExperiments extends Vue {
  readonly experimentContext = experimentModule.context(this.$store);

  headers = [
    {
      text: "Name",
      sortable: true,
      value: "name",
      align: "left",
    },
    {
      text: "Description",
      sortable: true,
      value: "description",
      align: "left",
    },
    {
      text: "Location",
      sortable: true,
      value: "location",
      align: "left",
    },
    {
      text: "Actions",
      value: "action",
      sortable: false,
    },
  ];

  get experiments() {
    return this.experimentContext.getters.experiments;
  }

  async mounted() {
    await this.experimentContext.actions.getExperiments();
  }

  async deleteExperiment(event, id: number) {
    if (self.confirm("Are you sure you want to delete the experiment?")) {
      await this.experimentContext.actions.deleteExperiment(id);
    }
  }
}
</script>
