<template>
  <div>
    <v-toolbar dense light>
      <v-toolbar-title>
        Manage Workflows
      </v-toolbar-title>
      <v-spacer></v-spacer>
      <v-btn small color="primary" to="/main/admin/workflows/create">Create Workflow</v-btn>
    </v-toolbar>

    <v-data-table
      :headers="headers"
      :items="workflows"
    >
      <template v-slot:item.action="{ item }">
        <v-tooltip bottom>
          <template v-slot:activator="{ on }">
            <v-btn v-on="on" icon :to="{name: 'main-admin-workflows-edit', params: {id: item.id}}">
              <v-icon>mdi-pencil</v-icon>
            </v-btn>
          </template>
          <span>Edit</span>
        </v-tooltip>
        <v-tooltip bottom>
          <template v-slot:activator="{ on }">
            <v-btn v-on="on" icon @click="deleteWorkflow($event, item.id)">
              <v-icon>mdi-delete</v-icon>
            </v-btn>
          </template>
          <span>Delete</span>
        </v-tooltip>
        <UploadButton :id="item.id"/>
      </template>
    </v-data-table>
  </div>
</template>

<script lang="ts">
  import UploadButton from '@/components/UploadButton.vue';
  import { workflowModule } from '@/modules/workflows';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: { UploadButton },
  })
  export default class AdminWorkflows extends Vue {
    readonly workflowContext = workflowModule.context(this.$store);

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
        text: 'Actions',
        value: 'action',
        sortable: false,
      },
    ];

    get workflows() {
      return this.workflowContext.getters.workflows;
    }

    async mounted() {
      await this.workflowContext.actions.getWorkflows();
    }

    async deleteWorkflow(event, id: number) {
      if (self.confirm('Are you sure you want to delete the workflow?')) {
        await this.workflowContext.actions.deleteWorkflow(id);
      }
    }
  }
</script>
