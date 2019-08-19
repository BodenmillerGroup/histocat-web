<template>
  <v-container fluid>
    <v-card class="ma-4 pa-4">
      <v-card-title primary-title>
        <div class="headline primary--text">Edit Workflow</div>
      </v-card-title>
      <v-card-text>
        <template>
          <v-form
            v-model="valid"
            ref="form"
            lazy-validation
          >
            <v-text-field
              label="Name"
              v-model="name"
              v-validate="'required'"
              data-vv-name="name"
              :error-messages="errors.collect('name')"
              required
            ></v-text-field>
            <v-text-field
              label="Description"
              v-model="description"
            ></v-text-field>
          </v-form>
        </template>
      </v-card-text>
      <v-card-actions>
        <v-spacer></v-spacer>
        <v-btn @click="cancel">Cancel</v-btn>
        <v-btn @click="reset">Reset</v-btn>
        <v-btn
          @click="submit"
          :disabled="!valid"
        >
          Save
        </v-btn>
      </v-card-actions>
    </v-card>
  </v-container>
</template>

<script lang="ts">
  import { workflowModule } from '@/modules/workflows';
  import { IWorkflowUpdate } from '@/modules/workflows/models';
  import { Component, Vue } from 'vue-property-decorator';

  @Component
  export default class EditWorkflow extends Vue {
    readonly workflowContext = workflowModule.context(this.$store);

    valid = true;
    name: string = '';
    description: string = '';

    get workflow() {
      return this.workflowContext.getters.getWorkflow(+this.$router.currentRoute.params.id);
    }

    async mounted() {
      this.reset();
    }

    reset() {
      this.name = '';
      this.description = '';
      this.$validator.reset();
    }

    cancel() {
      this.$router.back();
    }

    async submit() {
      if (await this.$validator.validateAll()) {
        const data: IWorkflowUpdate = {
          name: this.name,
        };
        if (this.description) {
          data.description = this.description;
        }
        await this.workflowContext.actions.updateWorkflow({ id: this.workflow!.id, data: data });
        this.$router.push('/main/admin/workflows');
      }
    }
  }
</script>
