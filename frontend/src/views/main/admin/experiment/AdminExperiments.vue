<template>
  <div>
    <v-toolbar light>
      <v-toolbar-title>
        Manage Experiments
      </v-toolbar-title>
      <v-spacer></v-spacer>
      <v-btn color="primary" to="/main/admin/experiments/create">Create Experiment</v-btn>
    </v-toolbar>
    <v-data-table :headers="headers" :items="experiments">
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
          <UploadButton :id="props.item.id" />
        </td>
      </template>
    </v-data-table>
  </div>
</template>

<script lang="ts">
  import { Component, Vue } from 'vue-property-decorator';
  import { readAdminExperiments } from '@/modules/experiment/getters';
  import { dispatchGetExperiments } from '@/modules/experiment/actions';
  import UploadButton from '@/components/UploadButton.vue';

  @Component({
    components: { UploadButton },
  })
  export default class AdminExperiments extends Vue {
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
        value: 'id',
      },
    ];

    get experiments() {
      return readAdminExperiments(this.$store);
    }

    async mounted() {
      await dispatchGetExperiments(this.$store);
    }
  }
</script>
