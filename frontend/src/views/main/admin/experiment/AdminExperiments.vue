<template>
  <div>
    <v-toolbar dense light>
      <v-toolbar-title>
        Manage Experiments
      </v-toolbar-title>
      <v-spacer></v-spacer>
      <v-btn small color="primary" to="/main/admin/experiments/create">Create Experiment</v-btn>
    </v-toolbar>

    <v-data-table
      :headers="headers"
      :items="experiments"
    >
      <template v-slot:item.action="{ item }">
        <v-tooltip bottom>
          <template v-slot:activator="{ on }">
            <v-btn v-on="on" icon :to="{name: 'main-admin-experiments-edit', params: {id: item.id}}">
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
        <UploadButton :id="item.id"/>
      </template>

      <template slot="items" slot-scope="props">
        <td>{{ props.item.name }}</td>
        <td>{{ props.item.description }}</td>
        <td>{{ props.item.location }}</td>
        <td class="justify-center layout px-0">
          <v-tooltip top>
            <span>Edit</span>
            <v-btn slot="activator" flat :to="{name: 'main-admin-experiments-edit', params: {id: props.item.id}}">
              <v-icon>mdi-pencil</v-icon>
            </v-btn>
          </v-tooltip>
          <v-tooltip top>
            <span>Delete</span>
            <v-btn slot="activator" flat @click="deleteExperiment($event, props.item.id)">
              <v-icon>mdi-delete</v-icon>
            </v-btn>
          </v-tooltip>
          <UploadButton :id="props.item.id"/>
        </td>
      </template>
    </v-data-table>
  </div>
</template>

<script lang="ts">
  import UploadButton from '@/components/UploadButton.vue';
  import { experimentModule } from '@/modules/experiment';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: { UploadButton },
  })
  export default class AdminExperiments extends Vue {
    readonly experimentContext = experimentModule.context(this.$store);

    headers = [
      {
        text: 'Name',
        sortable: true,
        value: 'name',
        align: 'left',
      },
      {
        text: 'Description',
        sortable: true,
        value: 'description',
        align: 'left',
      },
      {
        text: 'Location',
        sortable: true,
        value: 'location',
        align: 'left',
      },
      {
        text: 'Actions',
        value: 'action',
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
      const res = await this.$confirm('Do you really want to delete experiment?', { title: 'Warning' });
      if (res) {
        await this.experimentContext.actions.deleteExperiment(id);
      }
    }
  }
</script>
